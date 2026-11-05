# Optional linear pixel readout

`odd-digi-geometric-config.json` includes **disabled** `charge-quantization`
blocks for pixel volumes 16, 17 and 18. Native digitization remains the default;
all short- and long-strip settings are unchanged. This option requires ACTS with
support for `geometric.charge-quantization` and **`doMerge=true`**. Older ACTS
readers may silently ignore unknown JSON fields: check a JSON read/write round
trip in the runtime before producing data.

Select the configuration for the desired output:

- [Native](odd-digi-geometric-config.json): quantization disabled (default).
- [Four-bit P99.9](odd-digi-geometric-config-pixel-4bit-p999.json): quantization
  enabled for pixel volumes 16, 17 and 18. Only the three `enabled` flags differ
  from the native configuration; all strip settings remain unchanged.

ColliderML-Production accepts either file through its existing option:

```bash
--digi-config /path/to/ODD/config/odd-digi-geometric-config-pixel-4bit-p999.json
```

To return to native output, select `odd-digi-geometric-config.json`. Do not
enable quantization for entries with `digital=true`.

ACTS keeps its existing per-contribution charge smearing and thresholding, sums
accepted native contributions at identical channel addresses, then encodes each
aggregated channel exactly once. Neighboring-channel clustering and truth links
keep their existing behavior. Code weights determine the final measurement
position and the existing weighted uncertainty calculation. There is no added
noise, decoding, or second suppression threshold.

## Response and studied settings

For accepted activation `a`, `M = 2**bits - 1`:

```text
code = clamp(1 + floor((a - offset) / step), 1, M)
```

Boundary ties go into the upper bin. Zero is reserved for absent/suppressed
channels. Every retained channel receives a positive code. Four-bit output uses
codes 1–15; eight-bit output uses 1–255. Supported bit depths are 1–16.

The existing cell activation/value field contains the code for enabled pixel
entries. Other detector entries retain native values: do not treat the entire
cell table as either millimetres or codes without consulting its geometry and
configuration. Writer field types are unchanged; this is finite-precision
information, not packed four-bit storage. Native and code values can be studied
in separate runs using identical inputs and random seeds.

The enabled configuration uses the following fixed P99.9 four-bit parameters.
For P99.5 or P99, change only `step` in its three pixel blocks; these alternatives
do not require additional bundled detector configurations.

| Historical preset | `offset` | `step` |
| --- | ---: | ---: |
| P99.9 (example default) | 0.01 | 0.011105528959000351 |
| P99.5 | 0.01 | 0.009409294367857142 |
| P99 | 0.01 | 0.009057082750000002 |

These parameters came from 13,285,319 pixel cells in ColliderML ttbar pu200 v1,
run 0, calibration events 0–24 (validation events 25–49), using ACTS
`7cba36b173dd73a14342cc65043a43179a4b1dbd`, ODD v4.0.4
`b0992c148305224899a16d21fcd56406408bd393`. Input EDM4hep SHA256:
`7656dca207dfc96bb37c67ac524a6966d1a3ad10ebaaa4a9d3c0d493a874996f`;
digitization JSON SHA256:
`61ebc9fde798b9dcffcbda3b07ae9bce9352a4f00a2c1d4354dae6441b886a8d`.
Activation is a geometric path-length signal proxy in ACTS length units (mm),
including existing noise and merged contributions. The parameters are **not an
RD53 electronics calibration**, and must be validated against the chosen ODD
geometry and simulation inputs before recommending production use.

## Bit depth and overflow

Changing only `bits` preserves the bin width and increases the dynamic range.
For the same maximum-code entry point at a different bit depth, use:

```text
step_new = step_4bit * 14 / (2**bits_new - 2)
```

For example, use `bits: 8`, `offset: 0.01`, and
`step: 0.0006121157693937201` for the P99.9 range. This relation applies for
`bits_new >= 2`; one-bit positive readout has only code 1.

The maximum code begins at `offset + (M-1)*step`. Its nominal bin extends to
`offset + M*step`; only values beyond this latter boundary are the clipped
overflow tail. For P99.9 at four bits those points are approximately
0.1654774054 and 0.1765829344 mm. Therefore the fraction assigned code 15 and
the fraction beyond its nominal bin are different statistics.
