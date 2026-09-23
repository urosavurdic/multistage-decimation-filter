# Multistage decimation filter design

Designing a CIC decimation filter for R = 32 by exhaustive search: every
ordering of the prime factors, against every allocation of comb sections,
keeping the cheapest design that stays inside the passband deviation budget.

## The idea

Decimating by 32 in a single step needs a long, expensive filter. Doing it as a
cascade — one stage per prime factor — is much cheaper, because each stage runs
at a lower sample rate than the one before it. Two choices then drive the cost:

1. **The order the factors are applied in.** Decimating early means later
   stages run slower and cost less, but leaves less room to suppress aliasing.
2. **How many comb sections each stage gets.** More sections mean better
   attenuation and more arithmetic.

The search evaluates the operation count `A = K₁·R + Σ|Kₘ − Kₘ₊₁|·∏Rₖ + K_M`
for each candidate, sorts cheapest-first, and stops scanning as soon as the
cost exceeds the best design already found.

The frequency check only evaluates the bands that actually fold into the
passband under decimation — the response is unconstrained everywhere else, so
checking a uniform grid would waste most of the work.

For R = 32 every prime factor is 2, so there is exactly one ordering; the
ordering search only does real work when R has distinct prime factors.

## Running it

Open `decimation_filter_design.m` in MATLAB and run it. Needs the Signal
Processing Toolbox for `freqz`. It prints the chosen ordering, the comb-section
allocation, the operation count and the worst passband deviation.

| Path | What it is |
|---|---|
| `decimation_filter_design.m` | the search |
| `report.pdf` | written analysis |

## Credits

Coursework project with **Maja**, for Signal Processing 2.

## License

MIT — see [LICENSE](LICENSE).
