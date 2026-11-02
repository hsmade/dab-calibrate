#!/usr/bin/env python3
"""DAB Mode I reference-clock calibration from raw, unresampled IQ."""
import argparse
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
from numpy import fft

FS = 2048000
FRAME = 196608
N = 2048
CP = 504
SYMBOL = N + CP
CHANNELS = {'12C': 227360000, '11C': 220352000, '9C': 206352000}


def reference():
    """ETSI EN 300 401, Mode I phase-reference tables, in FFT-bin order.

    Numerical standard tables independently encoded here. Cross-reference:
    ODR-DabMod/src/PhaseReference.cpp (positive carriers, then negative).
    """
    h = np.array([
        [0,2,0,0,0,0,1,1,2,0,0,0,2,2,1,1],
        [0,3,2,3,0,1,3,0,2,1,2,3,2,3,3,0],
        [0,0,0,2,0,2,1,3,2,2,0,2,2,0,1,3],
        [0,1,2,1,0,3,3,2,2,3,2,1,2,1,3,2]])
    offsets = [3,1,1,1,2,2,1,0,2,2,3,3,0,2,1,3,3,3,3,0,3,0,1,1,
               1,2,0,1,3,2,2,3,2,1,2,3,1,2,3,3,2,2,2,1,1,3,1,2]
    rows = [0,3,2,1]*6 + [0,1,2,3]*6
    phases = np.concatenate([(np.tile(h[r], 2) + a) % 4
                             for r, a in zip(rows, offsets)])
    z = np.zeros(N, dtype=np.complex64)
    values = np.array([1,1j,-1,-1j], dtype=np.complex64)[phases]
    z[1:769], z[-768:] = values[:768], values[768:]
    return z


REF = reference()


class IQFile:
    def __init__(self, path, fmt):
        dtype, self.width = {'u8': ('u1', 2), 's16': ('<i2', 2),
                            'cf32': ('<c8', 1)}[fmt]
        self.fmt = fmt
        size = Path(path).stat().st_size
        unit = np.dtype(dtype).itemsize * self.width
        if not size or size % unit:
            raise ValueError('IQ file is empty or has an incomplete complex sample')
        self.data = np.memmap(path, dtype=dtype, mode='r')
        self.size = len(self.data) // self.width

    def read(self, lo, hi):
        raw = self.data[lo*self.width:hi*self.width]
        if self.fmt == 'cf32':
            x = np.array(raw, dtype=np.complex64)
            clipped = 0
        else:
            scale, midpoint = (128., 127.5) if self.fmt == 'u8' else (32768., 0.)
            clipped = np.count_nonzero((raw == (0 if self.fmt == 'u8' else -32768)) |
                                       (raw == (255 if self.fmt == 'u8' else 32767)))
            a = (raw.astype(np.float32).reshape(-1, 2) - midpoint) / scale
            x = a[:,0] + 1j*a[:,1]
        if not np.all(np.isfinite(x)):
            raise ValueError('IQ contains NaN or infinity')
        return x, clipped, len(raw)


def null_ends(x):
    """Power-envelope acquisition only; PRS subsequently confirms every candidate."""
    block = 64
    if len(x) < FRAME//2:
        return []
    p = np.mean(np.abs(x[:len(x)//block*block].reshape(-1, block))**2, axis=1)
    p = np.convolve(p, np.ones(5)/5, mode='same')
    level = np.median(p)
    if level <= 1e-12:
        return []
    # Noise fills the null on weak signals. Try progressively shallower dips;
    # cyclic-prefix coherence and the known PRS still validate each candidate.
    candidates = set()
    for fraction in (.30, .55, .75, .85):
        mask = p < level * fraction
        edges = np.diff(np.r_[False, mask, False].astype(int))
        starts, ends = np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)
        candidates.update(int(e*block) for s,e in zip(starts,ends)
                          if 1100 < (e-s)*block < 4200)
    return sorted(candidates)


def cp_offset(x, cp_start, count=1):
    # Trim the CP edges for timing uncertainty, sample-clock drift and multipath.
    starts = np.rint(cp_start + np.arange(count)*SYMBOL).astype(int)
    starts = starts[(starts >= 0) & (starts + N + 400 < len(x))]
    if not len(starts):
        raise ValueError('No complete OFDM symbols')
    idx = starts[:,None] + np.arange(120, 380)
    a, b = x[idx], x[idx+N]
    c = np.sum(b*np.conj(a), dtype=np.complex128)
    norm = np.sqrt(np.sum(abs(a)**2)*np.sum(abs(b)**2))
    return float(np.angle(c)*FS/(2*np.pi*N)), float(abs(c)/max(norm,1e-30))


def prs_measure(x, end, max_bins, previous=None, expected=None):
    cp = end
    s = cp + 240  # FFT window deliberately inside cyclic prefix.
    if s < 0 or s+N > len(x):
        return None
    fine, coherence = cp_offset(x, cp)
    if coherence < .12:
        return None
    y = fft.fft(x[s:s+N] * np.exp(-2j*np.pi*fine*np.arange(N)/FS))
    def search(shifts):
        shifted = y[(np.arange(N)[None,:] + shifts[:,None]) % N]
        cir = fft.ifft(shifted * np.conj(REF)[None,:], axis=1)
        power = abs(cir)**2
        # A correct PRS has concentrated channel impulse-response energy.
        scores = power.max(axis=1) / np.maximum(power.sum(axis=1), 1e-30)
        best = int(np.argmax(scores))
        return int(shifts[best]), float(scores[best]), power[best]
    shifts = np.arange(-max_bins, max_bins+1)
    if previous is not None:
        local = np.arange(max(-max_bins, previous-1), min(max_bins, previous+1)+1)
        k, score, p = search(local)
        if score < .08:
            k, score, p = search(shifts)
    else:
        k, score, p = search(shifts)
    if score < .04:
        return None
    peak = int(np.argmax(p))
    if expected is not None:
        # Keep one SFN path as relative strengths change. A path must remain
        # a local peak with at least 15% of the strongest path's power.
        peaks = np.flatnonzero((p > np.roll(p,1)) & (p >= np.roll(p,-1))
                               & (p > .15*p.max()))
        lags = np.where(peaks < N//2, peaks, peaks-N)
        distances = abs(s+lags-expected)
        if len(peaks) and distances.min() < 16:
            peak = int(peaks[np.argmin(distances)])
    # Signed circular lag; reject candidates outside plausible CP acquisition range.
    lag = peak if peak < N//2 else peak-N
    if not -100 < lag < 800:
        return None
    left, mid, right = p[(peak-1)%N], p[peak], p[(peak+1)%N]
    denominator = left - 2*mid + right
    sub = float(.5*(left-right)/denominator) if denominator < 0 else 0.
    useful = s + lag + np.clip(sub, -.5, .5)
    cp_actual = useful - CP
    fine2, coherence2 = cp_offset(x, cp_actual, 76)
    # Select the 1-kHz branch nearest the PRS solution (handles +/-500 wrap).
    total = k*1000 + fine
    total = fine2 + 1000*round((total-fine2)/1000)
    return {'position': float(useful), 'cfo_hz': float(total),
            'prs_score': score, 'cp_coherence': coherence2, 'coarse_bin': k}


def timing_fit(positions):
    positions = np.asarray(positions, dtype=float)
    if len(positions) < 8:
        raise ValueError('Need at least 8 confirmed DAB frames; use a stronger signal or longer capture')
    gaps = np.diff(positions)
    steps = np.rint(gaps / FRAME).astype(int)
    if np.any(steps < 1):
        raise ValueError('Ambiguous or duplicate frame detections')
    frame_numbers = np.r_[0, np.cumsum(steps)]
    keep = np.ones(len(positions), dtype=bool)
    for _ in range(5):
        b, a = np.polyfit(frame_numbers[keep], positions[keep]-positions[0], 1)
        residual = positions-positions[0] - (a+b*frame_numbers)
        med = np.median(residual[keep])
        mad = 1.4826*np.median(abs(residual[keep]-med))
        new = abs(residual-med) < max(2., 5*mad)
        if new.sum() < 8:
            raise ValueError('Timing fit rejected too many frames')
        if np.array_equal(new, keep):
            break
        keep = new
    b, a = np.polyfit(frame_numbers[keep], positions[keep]-positions[0], 1)
    residual = positions-positions[0] - (a+b*frame_numbers)
    rms = float(np.sqrt(np.mean(residual[keep]**2)))
    error = rms / np.sqrt(np.sum((frame_numbers[keep]-np.mean(frame_numbers[keep]))**2))
    # A discontinuity can be a USB loss or a switch to a different SFN path.
    jumps = int(np.count_nonzero(abs(np.diff(residual)) > 8))
    return {'sample_ppm': float((b/FRAME-1)*1e6),
            'fit_standard_error_ppm': float(error/FRAME*1e6),
            'timing_rms_samples': rms, 'samples_per_frame': float(b),
            'frames_used': int(keep.sum()), 'timing_jumps': jumps,
            'missing_frames': int(frame_numbers[-1]+1-len(positions)),
            'span_seconds': float((positions[-1]-positions[0])/FS)}


def analyze(path, fmt='u8', channel='12C', max_ppm=100., progress=True):
    iq = IQFile(path, fmt)
    if iq.size < FRAME*10:
        raise ValueError('Capture must contain at least one second of IQ')
    frames, clip, scalars = [], 0, 0
    chunk = FS*2
    overlap = FRAME*2
    last_position = -FRAME
    previous = None
    previous_ppm = 0.
    max_bins = int(math.ceil(CHANNELS[channel]*max_ppm/1e9)) + 2
    for lo in range(0, iq.size, chunk):
        if iq.size-lo < N+CP:
            break
        x, c, n = iq.read(lo, min(iq.size, lo+chunk+overlap))
        clip += c
        scalars += n
        x -= x.mean()
        for end in null_ends(x):
            if lo+end < last_position + FRAME*.5:
                continue
            expected = None
            if frames:
                gap = max(1, round((lo+end+CP-last_position)/FRAME))
                expected = last_position + gap*FRAME*(1+previous_ppm/1e6) - lo
            result = prs_measure(x, end, max_bins, previous, expected)
            if result is None:
                continue
            result['position'] += lo
            frames.append(result)
            last_position = result['position']
            previous = result['coarse_bin']
            previous_ppm = -result['cfo_hz']/CHANNELS[channel]*1e6
        if progress:
            print(f'\rAnalyzed {min(iq.size,lo+chunk)/FS:.1f}s; {len(frames)} PRS frames',
                  end='', file=sys.stderr, flush=True)
    if progress:
        print(file=sys.stderr)
    result = timing_fit([f['position'] for f in frames])
    cfo = np.array([f['cfo_hz'] for f in frames])
    median_cfo = float(np.median(cfo))
    rf_ppm = -median_cfo / CHANNELS[channel]*1e6
    warnings = []
    if result['span_seconds'] < 20:
        warnings.append('Short capture: use 30-60 seconds for calibration')
    if clip/max(scalars,1) > .001:
        warnings.append('More than 0.1% of IQ components reach ADC rails; reduce gain')
    if result['timing_jumps']:
        warnings.append('Timing discontinuities: possible sample loss or SFN path changes; recapture')
    if result['timing_rms_samples'] > 2:
        warnings.append('Unstable timing; clock estimate may be biased')
    if abs(rf_ppm-result['sample_ppm']) > 1:
        warnings.append('RF and sample estimates differ by >1 ppm; check IQ orientation, transmitter, and clock architecture')
    if abs(result['sample_ppm']) > max_ppm:
        warnings.append('Sample-clock estimate exceeds requested search range')
    if np.median([f['cp_coherence'] for f in frames]) < .15:
        warnings.append('Low cyclic-prefix coherence; RF estimate may be unreliable')
    if any(abs(f['coarse_bin']) >= max_bins for f in frames):
        warnings.append('Carrier estimate reaches search boundary; increase --max-ppm')
    result.update(channel=channel, frequency_hz=CHANNELS[channel], sample_rate_hz=FS,
                  input_format=fmt, frames_detected=len(frames),
                  rf_ppm=float(rf_ppm), carrier_offset_hz=median_cfo,
                  carrier_scatter_hz=float(1.4826*np.median(abs(cfo-median_cfo))),
                  median_prs_score=float(np.median([f['prs_score'] for f in frames])),
                  median_cp_coherence=float(np.median([f['cp_coherence'] for f in frames])),
                  clipping_fraction=float(clip/max(scalars,1)), warnings=warnings,
                  correction_recommended=not bool(result['timing_jumps'] or result['timing_rms_samples'] > 2
                                                  or abs(rf_ppm-result['sample_ppm']) > 1
                                                  or clip/max(scalars,1) > .001
                                                  or any(abs(f['coarse_bin']) >= max_bins for f in frames)
                                                  or np.median([f['cp_coherence'] for f in frames]) < .15
                                                  or abs(result['sample_ppm']) > max_ppm))
    return result


def capture(args):
    if not shutil.which('rtl_sdr'):
        raise ValueError('rtl_sdr executable not found; install the rtl-sdr package')
    path = Path(args.output)
    if path.exists() or path.with_suffix(path.suffix+'.json').exists():
        raise ValueError('Capture output or metadata already exists; choose a new path')
    path.parent.mkdir(parents=True, exist_ok=True)
    # Reserve output exclusively; rtl_sdr writes to our open descriptor.
    with path.open('xb') as out:
        # Odd sample count avoids legacy rtl_sdr exact-block termination bug.
        cmd = ['rtl_sdr','-d',args.device,'-f',str(CHANNELS[args.channel]),
               '-s',str(FS),'-p','0','-g',str(args.gain),
               '-b','262144','-n',str(round(args.seconds*FS) | 1),'-']
        print('Capturing with PPM correction = 0. Warm up the device first.', file=sys.stderr)
        try:
            process = subprocess.run(cmd, stdout=out, check=False, timeout=args.seconds+30)
        except subprocess.TimeoutExpired:
            raise ValueError(f'Capture exceeded time limit; partial IQ retained at {path}') from None
    if process.returncode:
        raise ValueError(f'rtl_sdr failed ({process.returncode}); partial capture retained at {path}')
    path.with_suffix(path.suffix+'.json').write_text(json.dumps(
        {'channel':args.channel,'sample_rate_hz':FS,'format':'u8','ppm_correction':0,
         'command':cmd}, indent=2)+'\n')
    return analyze(path, channel=args.channel, max_ppm=args.max_ppm, progress=not args.quiet)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    for name in ('analyze','capture'):
        p = sub.add_parser(name)
        p.add_argument('--channel', choices=CHANNELS, required=True)
        p.add_argument('--max-ppm', type=float, default=100)
        p.add_argument('--json', action='store_true', help='print machine-readable result')
        p.add_argument('--quiet', action='store_true')
        if name == 'analyze':
            p.add_argument('input', type=Path)
            p.add_argument('--format', choices=['u8','s16','cf32'], default='u8')
        else:
            p.add_argument('output', type=Path)
            p.add_argument('--seconds', type=float, default=60)
            p.add_argument('--gain', type=float, default=20)
            p.add_argument('--device', default='0')
    args = parser.parse_args()
    try:
        if not 0 < args.max_ppm <= 500:
            raise ValueError('--max-ppm must be between 0 and 500')
        if args.command == 'capture' and not 1 <= args.seconds <= 600:
            raise ValueError('--seconds must be between 1 and 600')
        result = capture(args) if args.command == 'capture' else analyze(
            args.input, args.format, args.channel, args.max_ppm, not args.quiet)
        if args.json:
            print(json.dumps(result, indent=2, allow_nan=False))
        else:
            print(f"{result['channel']} ({result['frequency_hz']/1e6:.3f} MHz)")
            print(f"Sample clock: {result['sample_ppm']:+.4f} ppm")
            print(f"RF estimate:  {result['rf_ppm']:+.4f} ppm")
            print(f"Carrier offset: {result['carrier_offset_hz']:+.2f} Hz")
            print(f"Frames: {result['frames_used']}/{result['frames_detected']}; "
                  f"span {result['span_seconds']:.1f}s; timing RMS {result['timing_rms_samples']:.3f} samples")
            print(f"Fit standard error: {result['fit_standard_error_ppm']:.4f} ppm (excludes systematic errors)")
            if result['correction_recommended']:
                print(f"Candidate librtlsdr correction: {round(result['sample_ppm']):+d} ppm; recapture to verify")
            else:
                print('No correction recommended: investigate quality warnings')
            for warning in result['warnings']:
                print(f'Warning: {warning}')
    except (OSError, ValueError) as error:
        print(f'Error: {error}', file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print('\nInterrupted; any partial capture is retained.', file=sys.stderr)
        return 130
    return 0


if __name__ == '__main__':
    sys.exit(main())
