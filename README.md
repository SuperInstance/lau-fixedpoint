# lau-fixedpoint

**Fixed-point arithmetic library for deterministic computation (Q16.16)**

A no-std-friendly, zero-dependency Rust library that wraps `i64` in a Q16.16 fixed-point format — 16 integer bits, 16 fractional bits — providing bit-identical results across every platform. No floating-point, no surprises.

---

## What This Does

Replaces `f64` arithmetic with a deterministic `FP` type and a companion `FPVec2` 2D vector. Includes trigonometric functions (`sin`, `cos`, `atan2`), square root (`fp_sqrt`), interpolation (`fp_lerp`), clamping, and a full set of operator overloads — all in pure integer math.

Use it when you need:
- **Replay determinism** — identical inputs always produce identical bit patterns, on every OS and CPU
- **Embedded / no-std targets** — no hardware FPU required
- **Networked game state** — fixed-point avoids cross-platform floating-point drift
- **Reproducible simulations** — physics, AI, or procedural generation that must replay identically

---

## Key Idea

Floating-point arithmetic is not associative and varies across hardware. Q16.16 fixed-point stores every number as `i64` with an implicit binary point between bit 16 and bit 15. Addition and subtraction are just integer ops; multiplication and division shift to compensate. The result is **bit-identical determinism** — the same calculation on any platform produces the same raw `i64`.

---

## Install

```toml
[dependencies]
lau-fixedpoint = "0.1"
```

Or:

```sh
cargo add lau-fixedpoint
```

Zero runtime dependencies. Rust 2021 edition.

---

## Quick Start

```rust
use lau_fixedpoint::*;

// Create fixed-point numbers
let a = FP::from_i64(3);       // 3.0
let b = FP::from_f64(1.5);     // 1.5
let c = FP::from_parts(2, 32768); // 2.5

// Arithmetic
assert_eq!((a + b).to_f64(), 4.5);
assert_eq!((a * b).to_f64(), 4.5);
assert_eq!((a / b).to_f64(), 2.0);

// Trigonometry
let sin_val = fp_sin(FP_HALF_PI);  // sin(π/2) ≈ 1.0
let cos_val = fp_cos(FP::ZERO);    // cos(0) ≈ 1.0
let angle = fp_atan2(FP::from_i64(1), FP::from_i64(1)); // ≈ π/4

// Square root
let root = fp_sqrt(FP::from_i64(144)); // ≈ 12.0

// Interpolation and clamping
let mid = fp_lerp(FP::from_i64(0), FP::from_i64(100), FP::HALF); // 50.0
let clamped = fp_clamp(FP::from_i64(150), FP::from_i64(0), FP::from_i64(100)); // 100

// 2D vectors
let v = FPVec2::new(FP::from_i64(3), FP::from_i64(4));
assert!((v.length().to_f64() - 5.0).abs() < 0.1);
let n = v.normalize();
assert!((n.length().to_f64() - 1.0).abs() < 0.05);

// Access raw parts
let val = FP::from_f64(3.5);
assert_eq!(val.integer_part(), 3);
assert_eq!(val.fraction_part(), 32768); // 0.5 in Q16.16
assert_eq!(val.raw(), 229376);          // 3.5 × 65536
```

---

## API Reference

### `FP` — Q16.16 Fixed-Point Number

A newtype wrapping `i64`. Derives `Clone`, `Copy`, `Debug`, `Default`, `Hash`, `Eq`, `PartialEq`, `Ord`, `PartialOrd`.

#### Constants

| Constant | Value | Description |
|---|---|---|
| `FP::ZERO` | 0 | Zero |
| `FP::ONE` | 65536 | 1.0 |
| `FP::HALF` | 32768 | 0.5 |
| `FP::MIN` | `i64::MIN` | Minimum representable |
| `FP::MAX` | `i64::MAX` | Maximum representable |
| `FP_PI` | ≈ 205887 | π |
| `FP_TWO_PI` | ≈ 411775 | 2π |
| `FP_HALF_PI` | ≈ 102944 | π/2 |
| `FP_E` | ≈ 178145 | e |

#### Constructors

| Method | Signature | Description |
|---|---|---|
| `from_i64` | `(i64) → FP` | Integer → fixed-point (shift left 16) |
| `from_f64` | `(f64) → FP` | Float → fixed-point (truncate) |
| `from_parts` | `(i32, u16) → FP` | Separate integer + fractional parts |
| `from_raw` | `(i64) → FP` | Wrap raw internal value directly |

#### Converters

| Method | Signature | Description |
|---|---|---|
| `to_i64` | `() → i64` | Truncate to integer (toward zero) |
| `to_f64` | `() → f64` | Convert to float (for display/debug) |
| `raw` | `() → i64` | Raw internal value |
| `integer_part` | `() → i32` | Integer component |
| `fraction_part` | `() → u16` | Raw 16-bit fractional component |

#### Operations

| Method | Signature | Description |
|---|---|---|
| `abs` | `() → FP` | Absolute value (saturating at MIN) |
| `scale` | `(i64) → FP` | Multiply by plain integer |

Operators: `+`, `-`, `*`, `/`, unary `-`. Uses saturating arithmetic for add/sub, widens to `i128` for mul/div to avoid overflow.

### Free Functions

| Function | Signature | Description |
|---|---|---|
| `fp_sqrt` | `(FP) → FP` | Newton's method, 8 iterations |
| `fp_sin` | `(FP) → FP` | Taylor series (5 terms), input in radians |
| `fp_cos` | `(FP) → FP` | `sin(x + π/2)` |
| `fp_atan2` | `(FP, FP) → FP` | CORDIC atan2, 16 iterations |
| `fp_lerp` | `(FP, FP, FP) → FP` | Linear interpolation: a + t(b − a) |
| `fp_clamp` | `(FP, FP, FP) → FP` | Clamp to [lo, hi] |
| `fp_min` | `(FP, FP) → FP` | Minimum |
| `fp_max` | `(FP, FP) → FP` | Maximum |

### `FPVec2` — 2D Fixed-Point Vector

```rust
pub struct FPVec2 { pub x: FP, pub y: FP }
```

| Method | Signature | Description |
|---|---|---|
| `new` | `(FP, FP) → FPVec2` | Constructor |
| `dot` | `(&FPVec2) → FP` | Dot product |
| `length` | `() → FP` | Magnitude via `fp_sqrt` |
| `normalize` | `() → FPVec2` | Unit vector (zero if length is zero) |
| `distance` | `(&FPVec2) → FP` | Euclidean distance |
| `lerp` | `(&FPVec2, FP) → FPVec2` | Component-wise interpolation |

Operators: `+`, `-`, unary `-`. Constant: `FPVec2::ZERO`.

---

## How It Works

### Q16.16 Format

Every `FP` value is an `i64` where the value equals `raw / 65536`. The top ~48 bits are the integer part (sign-extended), the bottom 16 bits are the fraction.

```
raw bits:  [sign | integer (47 bits) | fraction (16 bits)]
value:     sign × (integer_part + fraction_part / 65536)
```

### Multiplication

Two Q16.16 numbers multiplied give a Q32.32 result. The implementation widens to `i128`, multiplies, then right-shifts by 16 to get back to Q16.16:

```
result = ((a.raw as i128) × (b.raw as i128)) >> 16
```

### Division

Widens the dividend to `i128`, left-shifts by 16, then divides:

```
result = ((a.raw as i128) << 16) / (b.raw as i128)
```

Division by zero panics (no `NaN` in fixed-point).

### Square Root — Newton's Method

Starting from `guess = raw >> 1`, iterates 8 times:

```
guess = (guess + x / guess) / 2
```

Converges quickly for values up to ~2³¹ (the Q16.16 integer range).

### Sine — Taylor Series (5 terms)

```
sin(x) ≈ x − x³/3! + x⁵/5! − x⁷/7! + x⁹/9!
```

Input is first reduced modulo 2π to the range [−π, π]. The Taylor coefficients are precomputed in Q16.16. Accurate to within ~0.05 for typical inputs.

### atan2 — CORDIC

16 iterations of the CORDIC (COordinate Rotation DIgital Computer) algorithm using precomputed `atan(2^(−i))` angles. Produces results accurate to within ~0.1 radians. The quadrant is corrected based on the signs of x and y.

### Determinism

All operations are pure integer arithmetic — no FPU, no rounding modes, no extended precision. The determinism test computes 13 different operations twice and asserts bit-identical raw results.

---

## The Math

### Q-Format Representation

A Qm.n number stores `value × 2^n` as an integer. For Q16.16:
- Range: approximately −2⁴⁷ to 2⁴⁷ (with fractional precision)
- Resolution: 1/65536 ≈ 0.0000153

### Newton's Square Root

Iterating `g_{k+1} = (g_k + x/g_k) / 2` has quadratic convergence: the number of correct bits roughly doubles each iteration. 8 iterations from any positive starting guess gives excellent accuracy.

### Taylor Series for sin(x)

The Maclaurin series `sin(x) = Σ (−1)^k × x^(2k+1) / (2k+1)!` converges for all x. After reduction to [−π, π], 5 terms give accuracy better than |x^11 / 11!| ≈ 0.00003.

### CORDIC atan2

CORDIC rotates a vector (x, y) toward the x-axis using successive micro-rotations of angle `atan(2^(−i))`. After n iterations, the accumulated angle equals `atan2(y, x)` to within the remaining error.

### Saturating Arithmetic

Addition and subtraction use `i64::saturating_add`/`saturating_sub`, which clamp to `i64::MAX` or `i64::MIN` on overflow rather than wrapping — matching game/simulation expectations.

---

## Testing

34 tests covering:

- Construction from integers, floats, and parts
- Arithmetic operators (add, sub, mul, div, neg)
- Fractional arithmetic precision
- Comparison and ordering
- Square root (small, zero, large)
- Trigonometric functions (sin, cos, atan2)
- Utility functions (lerp, clamp, min, max)
- 2D vector operations (dot, length, normalize, distance, lerp, add, sub)
- Bit-identical determinism verification
- Constant accuracy (π, e)
- Division-by-zero panic

Run them:

```sh
cargo test
```

---

## License

MIT
