"""End-to-end waveform tests, including coarse CFO and sample-clock error."""
import sys
from pathlib import Path
import tempfile
import unittest

import numpy as np
from scipy import signal

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import dab_ppm as dab


def waveform(frames=32, ppm=20., channel='12C', snr=20., multipath=False):
    rng = np.random.default_rng(2026)
    bins = np.r_[np.arange(1,769), np.arange(1280,2048)]
    frame = np.zeros(dab.FRAME, dtype=np.complex64)
    for s in range(76):
        spec = np.zeros(dab.N, dtype=np.complex64)
        spec[bins] = np.exp(.5j*np.pi*rng.integers(0,4,len(bins)))
        if s == 0:
            spec = dab.reference()
        body = np.fft.ifft(spec).astype(np.complex64)*dab.N/np.sqrt(1536)
        start = 2656+s*dab.SYMBOL
        frame[start:start+dab.SYMBOL] = np.r_[body[-dab.CP:],body]
    x = np.tile(frame, frames)
    if multipath:
        x[37:] += (.25+.2j)*x[:-37]
    # 20ppm is exactly a 50001/50000 sample-rate ratio. Negative is inverse.
    if ppm == 20:
        x = signal.resample_poly(x,50001,50000)
    elif ppm == -20:
        x = signal.resample_poly(x,49999,50000)
    elif ppm != 0:
        raise ValueError('Test generator supports 0 and +/-20 ppm')
    cfo = -dab.CHANNELS[channel]*ppm/1e6
    x *= np.exp(2j*np.pi*cfo*np.arange(len(x))/dab.FS).astype(np.complex64)
    noise = 10**(-snr/20)/np.sqrt(2)
    x += noise*(rng.standard_normal(len(x))+1j*rng.standard_normal(len(x)))
    return x


class EstimatorTests(unittest.TestCase):
    def estimate(self,x,fmt='cf32',channel='12C'):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'test.iq'
            if fmt == 'cf32':
                x.astype('<c8').tofile(path)
            else:
                a = np.column_stack((x.real,x.imag))*.16
                np.clip(np.rint(a*128+127.5),0,255).astype('u1').tofile(path)
            return dab.analyze(path,fmt,channel,progress=False)

    def test_positive_clock_and_integer_cfo(self):
        r = self.estimate(waveform(multipath=True))
        self.assertAlmostEqual(r['sample_ppm'],20,delta=.08)
        self.assertAlmostEqual(r['rf_ppm'],20,delta=.08)
        self.assertGreater(r['frames_used'],25)
        self.assertTrue(r['correction_recommended'])

    def test_negative_clock_u8_channel_9c(self):
        r = self.estimate(waveform(ppm=-20,channel='9C'),fmt='u8',channel='9C')
        self.assertAlmostEqual(r['sample_ppm'],-20,delta=.08)
        self.assertAlmostEqual(r['rf_ppm'],-20,delta=.08)

    def test_zero_clock_channel_11c(self):
        r = self.estimate(waveform(ppm=0,channel='11C'),channel='11C')
        self.assertAlmostEqual(r['sample_ppm'],0,delta=.08)
        self.assertAlmostEqual(r['rf_ppm'],0,delta=.08)

    def test_sample_loss_is_flagged(self):
        x = waveform()
        x = np.r_[x[:len(x)//2], x[len(x)//2+64:]]
        r = self.estimate(x)
        self.assertGreater(r['timing_jumps'],0)
        self.assertFalse(r['correction_recommended'])

    def test_sfn_path_strength_switch_does_not_move_clock(self):
        x = waveform(ppm=0)
        echo = np.zeros_like(x)
        echo[62:] = x[:-62]
        amplitude = np.full(len(x), .55, dtype=np.float32)
        amplitude[dab.FRAME*8:dab.FRAME*22] = 1.4
        x += echo*amplitude
        r = self.estimate(x)
        self.assertAlmostEqual(r['sample_ppm'],0,delta=.08)
        self.assertEqual(r['timing_jumps'],0)
        self.assertTrue(r['correction_recommended'])

    def test_weak_signal_with_noise_filled_null(self):
        r = self.estimate(waveform(snr=-3))
        self.assertAlmostEqual(r['sample_ppm'],20,delta=.15)
        self.assertAlmostEqual(r['rf_ppm'],20,delta=.15)
        self.assertGreater(r['frames_used'],20)
        self.assertTrue(r['correction_recommended'])

    def test_noise_does_not_lock(self):
        rng = np.random.default_rng(10)
        x = (rng.normal(size=dab.FS)+1j*rng.normal(size=dab.FS)).astype('c8')
        with self.assertRaises(ValueError):
            self.estimate(x)

    def test_missing_frames_preserve_frame_index(self):
        k = np.r_[np.arange(30),np.arange(33,80)]
        p = 1234. + k*dab.FRAME*1.000020
        r = dab.timing_fit(p)
        self.assertAlmostEqual(r['sample_ppm'],20,places=6)
        self.assertEqual(r['missing_frames'],3)

    def test_reference_layout_and_standard_values(self):
        z = dab.reference()
        self.assertEqual(np.count_nonzero(z),1536)
        self.assertEqual(z[0],0)
        self.assertEqual(z[1],-1j)
        self.assertEqual(z[-768],1j)
        self.assertEqual(z[33],1j)

    def test_truncated_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'bad.iq'
            path.write_bytes(b'123')
            with self.assertRaises(ValueError):
                dab.IQFile(path,'u8')


if __name__ == '__main__':
    unittest.main()
