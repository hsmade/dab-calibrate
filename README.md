# DAB PPM calibrator

Measure an RTL-SDR's sample-clock error and RF tuning error using DAB Mode I
on 12C (227.360 MHz), 11C (220.352 MHz), or 9C (206.352 MHz).
No audio decoder, GPS receiver, or computer wall-clock reference is needed.

## Docker

The image includes Python, NumPy, and a pinned RTL-SDR Blog driver with V4
support. It supports `linux/amd64` and `linux/arm64` (including 64-bit Raspberry Pi).

Build locally:

```sh
docker build -t dab-calibrate .
docker run --rm dab-calibrate --help
```

Capture on the machine with the USB receiver attached:

```sh
mkdir -p data
# Find the receiver's bus/device address using lsusb; substitute its device node.
docker run --rm \
  --device /dev/bus/usb/003/002 \
  -v "$PWD/data:/data" \
  dab-calibrate capture /data/12C.u8 --channel 12C --device 0 --seconds 60 --gain 29 --json
```

Choose `--device` using the enumeration printed by `rtl_sdr`; other receivers
may still appear in the list even when their USB nodes are not mapped. Stop any
other program using the selected receiver before capturing. Captures are written to the
mounted directory and may be owned by root; `--user "$(id -u):$(id -g)"` can be used
when that UID/GID has access to both the USB node and output directory. The image
runs as root by default to support USB access. No `--privileged` flag is required.

Analyze an existing recording, without USB access:

```sh
docker run --rm -v "$PWD/data:/data:ro" \
  dab-calibrate analyze /data/12C.u8 --channel 12C --json
```

After the publishing workflow has run, substitute
`ghcr.io/hsmade/dab-calibrate:latest` for `dab-calibrate` in these commands.
The same entrypoint supports both `capture` and `analyze`.

## GitHub Container Registry

`.github/workflows/docker.yml` runs waveform tests and builds/smoke-tests the
image on pull requests and pushes to `main`. Successful pushes to `main`, Git tag
pushes, and manual workflow runs publish both architectures to
`ghcr.io/<repository-owner>/<repository-name>` using `GITHUB_TOKEN`.

- Main or manual build: `:latest`.
- Git tag `v1.2.3`: both `:latest` and `:v1.2.3`.
- Pull requests: tests and local image build only, with no publishing credentials.

Use Docker-compatible Git tag names for exact matching. Docker metadata normalizes
characters that are invalid in image tags. A tagged build updates `latest` too,
including when the tag refers to an older commit. No extra registry secret is needed.
Actions are pinned to commit SHAs, and the RTL driver is pinned to a source commit
in `Dockerfile`. The driver license, corresponding source archive, and source commit are included in the image at
`/usr/local/share/doc/rtl-sdr/`.

## Run without Docker

Python 3.10+ and NumPy are required. Install into a virtual environment:

```sh
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

Warm up the dongle for 10–15 minutes. Use a stationary antenna and choose a
channel with strong reception. Capture and analyze device 0:

```sh
.venv/bin/python dab_ppm.py capture 12C.u8 --channel 12C --device 0 --seconds 60 --gain 29
```

`capture` invokes `rtl_sdr` from PATH, with PPM correction explicitly set to zero.
It refuses to overwrite existing output, retains IQ, and writes a capture metadata
sidecar. It uses an odd sample count to avoid a legacy rtl_sdr bug where a capture
whose length is an exact number of USB blocks runs indefinitely.

Alternatively, record elsewhere and analyze locally:

```sh
rtl_sdr -d 0 -f 227360000 -s 2048000 -p 0 -g 29 -n 122880001 12C.u8
.venv/bin/python dab_ppm.py analyze 12C.u8 --channel 12C
.venv/bin/python dab_ppm.py analyze 12C.u8 --channel 12C --json > result.json
```

For remote receivers, copy `dab_ppm.py` to the host and run `capture` there,
then transfer its small JSON result. Capture IQ to the host's local disk rather
than piping it over SSH: network backpressure can block USB processing and lose
samples. The production tool only needs NumPy on the remote host.

The raw file MUST contain **unresampled, 2.048 MS/s IQ**. Supported formats:

- `--format u8`: unsigned interleaved I,Q bytes (rtl_sdr output; default).
- `--format s16`: little-endian signed interleaved I,Q 16-bit samples.
- `--format cf32`: little-endian complex float32, I then Q.

No WAV, SigMF-header parsing, or automatic sample-rate conversion is performed.
The channel argument must match the recording. Default acquisition range is
+/-100 ppm; `--max-ppm 200` expands it. Use a V4-compatible driver for Blog V4
hardware. Older distribution drivers may open the device successfully yet produce
no usable signal. A private driver can be selected using PATH.

## Interpretation and verification

`Sample clock: +20 ppm` means the receiver clock is fast relative to DAB frame
timing. `RF estimate: +20 ppm` means the frequency error is consistent with a
fast tuner reference, assuming conventional IQ orientation. For a shared clock,
these should agree within measurement and transmitter errors.

With standard librtlsdr, the candidate correction for this example is `-p 20`
(positive 20, not negative). Apply it only after inspecting quality warnings.
Recapture to verify that residual sample and RF PPM approach zero. An already
corrected recording reports **residual** error; the displayed candidate assumes
the recording was made with `-p 0`. This tool does not change a persistent device
setting or EEPROM. Integer-PPM APIs cannot apply a fractional correction exactly.

## How it works

1. Detect the DAB null interval in the power envelope.
2. Estimate fractional carrier offset by averaging cyclic-prefix correlations
   across the available frame (modulo 1 kHz).
3. Search integer carrier shifts using the known ETSI phase-reference symbol.
   The correct shift produces a concentrated channel impulse response.
4. Interpolate an impulse-response peak for each frame's timing, retaining the
   same SFN path while its strength changes (when it remains detectable).
5. Fit sample position versus frame number, rejecting timing outliers and allowing
   missed detections. The nominal frame is 196,608 samples / 96 ms.
6. Refine carrier offset using cyclic prefixes across the available frame.

Formulas, for zero initial PPM correction:

```text
sample_ppm = 1e6 * (fitted_samples_per_frame / 196608 - 1)
rf_ppm     = -1e6 * measured_baseband_offset_hz / nominal_channel_hz
```

Raw input is memory-mapped and processed in overlapping two-second chunks, so
memory use does not scale with recording length. Results include correlation
quality, timing residuals, carrier scatter, clipping, and possible discontinuities.

## Limits

- Accuracy is relative to the transmitter clock and carrier; this is not a
  traceable absolute frequency standard. Test multiple ensembles for consistency.
- The fit standard error describes random timing residuals only. It excludes
  transmitter error, multipath, temperature drift, and systematic interpolation bias.
- A changing strongest SFN path can mimic timing drift. Abrupt timing jumps suppress
  a correction recommendation; slow changes can remain undetected.
- Lost USB samples can bias results. Large non-frame timing jumps are flagged;
  whole-frame loss and very small discontinuities cannot always be identified.
- Progressively shallower null thresholds support weak signals, but very low SNR
  or TII-filled nulls can still prevent acquisition. No confirmed PRS lock
  produces an error, rather than a PPM estimate from unverified power dips.
- Sample-rate and tuner clocks may be separate on modified hardware. RF PPM is
  approximate because nominal sample rate scales the frequency measurement.
- Clipping checks count ADC rail components; floating-point files have no rail check.
- The program assumes ordinary I+jQ. Conjugated IQ is not auto-detected.

## Tests

```sh
python3 -m pip install -r requirements-test.txt
python3 -m unittest discover -s tests -v
```

The waveform tests generate full Mode I frames with PRS, cyclic prefixes and
random OFDM payload, inject +/-20 ppm resampling, multi-kHz tuning errors, AWGN,
multipath, quantization, noise-filled nulls at -3 and -10 dB SNR, and dropped samples. Noise-only input must fail to lock.
These are physical-layer tests, not FIC or audio decoding tests.

## References

- ETSI EN 300 401, Mode I frame timing and phase-reference numerical tables:
  https://www.etsi.org/deliver/etsi_i_ets/300400_300499/300401/02_60/ets_300401e02p.pdf
- PRS numerical table cross-check:
  https://github.com/Opendigitalradio/ODR-DabMod/blob/master/src/PhaseReference.cpp
- PRS timing as a sample-clock measurement:
  https://qirx.softsyst.com/QIRXCalibrate?sequenceNo=2
- RTL-SDR Blog V4 driver requirements: https://www.rtl-sdr.com/v4/
