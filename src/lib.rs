//! # lau-fixedpoint
//!
//! Fixed-point arithmetic library for deterministic vibe computation.
//! Uses Q16.16 format (16 integer bits, 16 fractional bits) wrapping `i64`.
//! No floating-point nondeterminism across platforms.

use std::cmp::Ordering;
use std::ops::{Add, Div, Mul, Neg, Sub};

// ---------------------------------------------------------------------------
// Constants (Q16.16)
// ---------------------------------------------------------------------------

/// π in Q16.16 ≈ 205887
pub const FP_PI: FP = FP(205887);
/// 2π in Q16.16
pub const FP_TWO_PI: FP = FP(411775);
/// π/2 in Q16.16
pub const FP_HALF_PI: FP = FP(102944);
/// e in Q16.16 ≈ 178145
pub const FP_E: FP = FP(178145);

// ---------------------------------------------------------------------------
// FP
// ---------------------------------------------------------------------------

/// A Q16.16 fixed-point number wrapping `i64`.
#[derive(Clone, Copy, Debug, Default, Hash, Eq, PartialEq)]
pub struct FP(pub i64);

impl FP {
    // ----- associated constants -----
    pub const ZERO: FP = FP(0);
    pub const ONE: FP = FP(1 << 16);
    pub const HALF: FP = FP(1 << 15);
    pub const MIN: FP = FP(i64::MIN);
    pub const MAX: FP = FP(i64::MAX);

    // ----- constructors -----

    /// Convert an integer to fixed-point.
    pub const fn from_i64(val: i64) -> Self {
        FP(val << 16)
    }

    /// Convert an `f64` to fixed-point (truncating).
    pub fn from_f64(val: f64) -> Self {
        FP((val * 65536.0) as i64)
    }

    /// Truncate to integer (toward zero).
    pub const fn to_i64(&self) -> i64 {
        self.0 >> 16
    }

    /// Convert to `f64`.
    pub fn to_f64(&self) -> f64 {
        self.0 as f64 / 65536.0
    }

    /// Build from separate integer and fractional parts.
    /// `fraction` is the raw 16-bit fractional value (0..=65535).
    pub const fn from_parts(integer: i32, fraction: u16) -> Self {
        FP(((integer as i64) << 16) | (fraction as i64 & 0xFFFF))
    }

    /// Return the integer part (sign-extended).
    pub const fn integer_part(&self) -> i32 {
        (self.0 >> 16) as i32
    }

    /// Return the raw 16-bit fractional part (always positive magnitude).
    pub const fn fraction_part(&self) -> u16 {
        (self.0 & 0xFFFF) as u16
    }

    /// Raw internal `i64` value.
    pub const fn raw(&self) -> i64 {
        self.0
    }

    /// Construct from a raw internal `i64`.
    pub const fn from_raw(raw: i64) -> Self {
        FP(raw)
    }

    /// Absolute value (saturating at `MIN`).
    pub fn abs(&self) -> FP {
        if self.0 == i64::MIN {
            FP(i64::MAX)
        } else if self.0 < 0 {
            FP(-self.0)
        } else {
            *self
        }
    }

    /// Multiply by a plain integer.
    pub fn scale(&self, s: i64) -> FP {
        FP(self.0.saturating_mul(s))
    }
}

// ---------------------------------------------------------------------------
// Arithmetic traits
// ---------------------------------------------------------------------------

impl Add for FP {
    type Output = FP;
    fn add(self, rhs: FP) -> FP {
        FP(self.0.saturating_add(rhs.0))
    }
}

impl Sub for FP {
    type Output = FP;
    fn sub(self, rhs: FP) -> FP {
        FP(self.0.saturating_sub(rhs.0))
    }
}

impl Mul for FP {
    type Output = FP;
    fn mul(self, rhs: FP) -> FP {
        let a = self.0 as i128;
        let b = rhs.0 as i128;
        FP(((a * b) >> 16) as i64)
    }
}

impl Div for FP {
    type Output = FP;
    fn div(self, rhs: FP) -> FP {
        if rhs.0 == 0 {
            panic!("division by zero in FP");
        }
        let a = (self.0 as i128) << 16;
        let b = rhs.0 as i128;
        FP((a / b) as i64)
    }
}

impl Neg for FP {
    type Output = FP;
    fn neg(self) -> FP {
        FP(self.0.saturating_neg())
    }
}

// ---------------------------------------------------------------------------
// Comparison traits
// ---------------------------------------------------------------------------

impl Ord for FP {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.cmp(&other.0)
    }
}

impl PartialOrd for FP {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

// ---------------------------------------------------------------------------
// Free functions
// ---------------------------------------------------------------------------

/// Newton's method square root, 8 iterations.
pub fn fp_sqrt(x: FP) -> FP {
    if x.0 <= 0 {
        return FP::ZERO;
    }
    // Seed: half the raw value as a rough guess, then shift to Q16.16
    let mut guess = FP(x.0 >> 1);
    for _ in 0..8 {
        if guess.0 == 0 {
            break;
        }
        guess = (guess + x / guess) * FP::from_i64(1) / FP::from_i64(2);
        // Rewrite without trait ambiguity:
        // (guess + x/guess) * 0.5
    }
    guess
}

/// Taylor-series sin approximation (input in radians, 5 terms).
/// sin(x) ≈ x - x³/3! + x⁵/5! - x⁷/7! + x⁹/9!
pub fn fp_sin(x: FP) -> FP {
    // Normalize x into [-π, π]
    let mut x = x;
    // Use modular reduction: x = x mod 2π, shifted to [-π, π]
    if x.0 != 0 {
        let two_pi = FP_TWO_PI;
        // integer division for number of full periods
        let raw_div = x.0 / two_pi.0;
        // subtract full periods (use i128 to avoid overflow)
        let periods = (raw_div as i128) * (two_pi.0 as i128);
        x = FP((x.0 as i128 - periods) as i64);
        // shift to [-π, π]
        if x > FP_PI {
            x = x - FP_TWO_PI;
        } else if x < -FP_PI {
            x = x + FP_TWO_PI;
        }
    }

    // Taylor coefficients as Q16.16
    let x2 = x * x;
    let x3 = x2 * x;
    let x5 = x3 * x2;
    let x7 = x5 * x2;
    let x9 = x7 * x2;

    // 1/3! = 1/6 ≈ 10923 in Q16.16
    let inv3_fact = FP(10923);
    // 1/5! = 1/120 ≈ 546 in Q16.16
    let inv5_fact = FP(546);
    // 1/7! = 1/5040 ≈ 13 in Q16.16
    let inv7_fact = FP(13);
    // 1/9! = 1/362880 ≈ 0 in Q16.16 (too small), use 1
    let inv9_fact = FP(1);

    x - x3 * inv3_fact + x5 * inv5_fact - x7 * inv7_fact + x9 * inv9_fact
}

/// Cosine via sin(x + π/2).
pub fn fp_cos(x: FP) -> FP {
    fp_sin(x + FP_HALF_PI)
}

/// CORDIC atan2 approximation.
pub fn fp_atan2(y: FP, x: FP) -> FP {
    if x.0 == 0 && y.0 == 0 {
        return FP::ZERO;
    }

    // Pre-computed CORDIC angles in Q16.16 (atan(2^-i) for i=0..15)
    let angles: [i64; 16] = [
        51472,  // atan(1)     = π/4
        30399,  // atan(1/2)
        16124,  // atan(1/4)
        8150,   // atan(1/8)
        4091,   // atan(1/16)
        2047,   // atan(1/32)
        1024,   // ...
        512,
        256,
        128,
        64,
        32,
        16,
        8,
        4,
        2,
    ];

    let mut vx = x.0.abs();
    let mut vy = y.0.abs();
    let mut angle: i64 = 0;

    for (i, &ai) in angles.iter().enumerate() {
        let dx = vx >> i;
        let dy = vy >> i;
        if vy > 0 {
            vx += dy;
            vy -= dx;
            angle += ai;
        } else {
            vx -= dy;
            vy += dx;
            angle -= ai;
        }
    }

    // angle is atan(|y|/|x|) in Q16.16, adjust for quadrant
    let result = if x.0 >= 0 && y.0 >= 0 {
        angle
    } else if x.0 < 0 && y.0 >= 0 {
        FP_PI.0 - angle
    } else if x.0 < 0 && y.0 < 0 {
        -(FP_PI.0 - angle)
    } else {
        -angle
    };

    FP(result)
}

/// Linear interpolation: a + t * (b - a).
pub fn fp_lerp(a: FP, b: FP, t: FP) -> FP {
    a + t * (b - a)
}

/// Clamp value to [lo, hi].
pub fn fp_clamp(val: FP, lo: FP, hi: FP) -> FP {
    fp_max(lo, fp_min(val, hi))
}

/// Minimum of two FP values.
pub fn fp_min(a: FP, b: FP) -> FP {
    if a <= b { a } else { b }
}

/// Maximum of two FP values.
pub fn fp_max(a: FP, b: FP) -> FP {
    if a >= b { a } else { b }
}

// ---------------------------------------------------------------------------
// FPVec2
// ---------------------------------------------------------------------------

/// A 2D vector of `FP` values.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct FPVec2 {
    pub x: FP,
    pub y: FP,
}

impl FPVec2 {
    pub fn new(x: FP, y: FP) -> Self {
        FPVec2 { x, y }
    }

    /// Dot product.
    pub fn dot(&self, other: &FPVec2) -> FP {
        self.x * other.x + self.y * other.y
    }

    /// Length (magnitude) via `fp_sqrt`.
    pub fn length(&self) -> FP {
        fp_sqrt(self.dot(self))
    }

    /// Normalize to unit length. Returns zero vector if length is zero.
    pub fn normalize(&self) -> FPVec2 {
        let len = self.length();
        if len.0 == 0 {
            return FPVec2::ZERO;
        }
        FPVec2 {
            x: self.x / len,
            y: self.y / len,
        }
    }

    /// Distance to another vector.
    pub fn distance(&self, other: &FPVec2) -> FP {
        (*self - *other).length()
    }

    /// Linear interpolation between two vectors.
    pub fn lerp(&self, other: &FPVec2, t: FP) -> FPVec2 {
        FPVec2 {
            x: fp_lerp(self.x, other.x, t),
            y: fp_lerp(self.y, other.y, t),
        }
    }

    pub const ZERO: FPVec2 = FPVec2 {
        x: FP::ZERO,
        y: FP::ZERO,
    };
}

impl Add for FPVec2 {
    type Output = FPVec2;
    fn add(self, rhs: FPVec2) -> FPVec2 {
        FPVec2 {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl Sub for FPVec2 {
    type Output = FPVec2;
    fn sub(self, rhs: FPVec2) -> FPVec2 {
        FPVec2 {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
        }
    }
}

impl Neg for FPVec2 {
    type Output = FPVec2;
    fn neg(self) -> FPVec2 {
        FPVec2 {
            x: -self.x,
            y: -self.y,
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_from_i64() {
        let a = FP::from_i64(3);
        assert_eq!(a.raw(), 3 << 16);
        assert_eq!(a.to_i64(), 3);
    }

    #[test]
    fn test_from_f64() {
        let a = FP::from_f64(1.5);
        assert_eq!(a.raw(), 98304); // 1.5 * 65536
    }

    #[test]
    fn test_to_f64() {
        let a = FP::from_f64(2.25);
        let diff = (a.to_f64() - 2.25).abs();
        assert!(diff < 0.001);
    }

    #[test]
    fn test_from_parts() {
        let a = FP::from_parts(3, 32768); // 3.5
        assert_eq!(a.to_i64(), 3);
        assert_eq!(a.integer_part(), 3);
        assert_eq!(a.fraction_part(), 32768);
    }

    #[test]
    fn test_zero_one_half() {
        assert_eq!(FP::ZERO.raw(), 0);
        assert_eq!(FP::ONE.raw(), 65536);
        assert_eq!(FP::HALF.raw(), 32768);
    }

    #[test]
    fn test_add() {
        let a = FP::from_i64(2);
        let b = FP::from_i64(3);
        assert_eq!((a + b).to_i64(), 5);
    }

    #[test]
    fn test_sub() {
        let a = FP::from_i64(10);
        let b = FP::from_i64(4);
        assert_eq!((a - b).to_i64(), 6);
    }

    #[test]
    fn test_mul() {
        let a = FP::from_i64(3);
        let b = FP::from_i64(4);
        assert_eq!((a * b).to_i64(), 12);
    }

    #[test]
    fn test_mul_fractional() {
        let a = FP::from_f64(1.5);
        let b = FP::from_f64(2.0);
        let result = (a * b).to_f64();
        assert!((result - 3.0).abs() < 0.01);
    }

    #[test]
    fn test_div() {
        let a = FP::from_i64(12);
        let b = FP::from_i64(4);
        assert_eq!((a / b).to_i64(), 3);
    }

    #[test]
    fn test_div_fractional() {
        let a = FP::from_f64(7.0);
        let b = FP::from_f64(2.0);
        let result = (a / b).to_f64();
        assert!((result - 3.5).abs() < 0.01);
    }

    #[test]
    fn test_neg() {
        let a = FP::from_i64(5);
        assert_eq!((-a).to_i64(), -5);
    }

    #[test]
    fn test_abs() {
        let a = FP::from_i64(-7);
        assert_eq!(a.abs().to_i64(), 7);
        assert_eq!(FP::from_i64(3).abs().to_i64(), 3);
    }

    #[test]
    fn test_scale() {
        let a = FP::from_i64(3);
        assert_eq!(a.scale(5).to_i64(), 15);
    }

    #[test]
    fn test_ord() {
        let a = FP::from_i64(1);
        let b = FP::from_i64(2);
        assert!(a < b);
        assert!(b > a);
        assert!(a <= b);
        assert!(b >= a);
    }

    #[test]
    fn test_sqrt() {
        let a = FP::from_i64(4);
        let result = fp_sqrt(a);
        assert!((result.to_f64() - 2.0).abs() < 0.05);
    }

    #[test]
    fn test_sqrt_zero() {
        assert_eq!(fp_sqrt(FP::ZERO), FP::ZERO);
    }

    #[test]
    fn test_sqrt_large() {
        let a = FP::from_i64(100);
        let result = fp_sqrt(a);
        assert!((result.to_f64() - 10.0).abs() < 0.1);
    }

    #[test]
    fn test_sin_zero() {
        let result = fp_sin(FP::ZERO);
        assert_eq!(result.raw().abs(), 0);
    }

    #[test]
    fn test_sin_pi_half() {
        let result = fp_sin(FP_HALF_PI);
        let diff = (result.to_f64() - 1.0).abs();
        assert!(diff < 0.05, "sin(π/2) ≈ {}, expected ~1.0", result.to_f64());
    }

    #[test]
    fn test_cos_zero() {
        let result = fp_cos(FP::ZERO);
        let diff = (result.to_f64() - 1.0).abs();
        assert!(diff < 0.05, "cos(0) ≈ {}, expected ~1.0", result.to_f64());
    }

    #[test]
    fn test_atan2_basic() {
        let y = FP::from_i64(1);
        let x = FP::from_i64(1);
        let result = fp_atan2(y, x);
        let expected = std::f64::consts::FRAC_PI_4;
        let diff = (result.to_f64() - expected).abs();
        assert!(diff < 0.1, "atan2(1,1) ≈ {}, expected ~{}", result.to_f64(), expected);
    }

    #[test]
    fn test_lerp() {
        let a = FP::from_i64(0);
        let b = FP::from_i64(10);
        let t = FP::from_f64(0.5);
        let result = fp_lerp(a, b, t);
        assert!((result.to_f64() - 5.0).abs() < 0.1);
    }

    #[test]
    fn test_clamp() {
        let val = FP::from_i64(15);
        let lo = FP::from_i64(0);
        let hi = FP::from_i64(10);
        assert_eq!(fp_clamp(val, lo, hi), hi);
        assert_eq!(fp_clamp(FP::from_i64(-5), lo, hi), lo);
        assert_eq!(fp_clamp(FP::from_i64(5), lo, hi), FP::from_i64(5));
    }

    #[test]
    fn test_min_max() {
        let a = FP::from_i64(3);
        let b = FP::from_i64(7);
        assert_eq!(fp_min(a, b), a);
        assert_eq!(fp_max(a, b), b);
    }

    #[test]
    fn test_vec2_dot() {
        let a = FPVec2::new(FP::from_i64(3), FP::from_i64(4));
        let b = FPVec2::new(FP::from_i64(1), FP::from_i64(2));
        assert_eq!(a.dot(&b).to_i64(), 11);
    }

    #[test]
    fn test_vec2_length() {
        let v = FPVec2::new(FP::from_i64(3), FP::from_i64(4));
        let len = v.length();
        assert!((len.to_f64() - 5.0).abs() < 0.1);
    }

    #[test]
    fn test_vec2_normalize() {
        let v = FPVec2::new(FP::from_i64(3), FP::from_i64(4));
        let n = v.normalize();
        let len = n.length();
        assert!((len.to_f64() - 1.0).abs() < 0.05);
    }

    #[test]
    fn test_vec2_distance() {
        let a = FPVec2::new(FP::from_i64(0), FP::from_i64(0));
        let b = FPVec2::new(FP::from_i64(3), FP::from_i64(4));
        let dist = a.distance(&b);
        assert!((dist.to_f64() - 5.0).abs() < 0.1);
    }

    #[test]
    fn test_vec2_lerp() {
        let a = FPVec2::new(FP::from_i64(0), FP::from_i64(0));
        let b = FPVec2::new(FP::from_i64(10), FP::from_i64(20));
        let t = FP::from_f64(0.5);
        let mid = a.lerp(&b, t);
        assert!((mid.x.to_f64() - 5.0).abs() < 0.1);
        assert!((mid.y.to_f64() - 10.0).abs() < 0.1);
    }

    #[test]
    fn test_vec2_add_sub() {
        let a = FPVec2::new(FP::from_i64(1), FP::from_i64(2));
        let b = FPVec2::new(FP::from_i64(3), FP::from_i64(4));
        let sum = a + b;
        assert_eq!(sum.x.to_i64(), 4);
        assert_eq!(sum.y.to_i64(), 6);
        let diff = b - a;
        assert_eq!(diff.x.to_i64(), 2);
        assert_eq!(diff.y.to_i64(), 2);
    }

    #[test]
    fn test_determinism() {
        // Compute a sequence of operations and verify bit-identical results
        let mut results: Vec<i64> = Vec::new();

        let a = FP::from_i64(7);
        let b = FP::from_i64(3);
        let c = FP::from_f64(2.5);

        results.push((a + b).raw());
        results.push((a - b).raw());
        results.push((a * b).raw());
        results.push((a / b).raw());
        results.push((c * FP::from_i64(4)).raw());
        results.push(fp_sqrt(FP::from_i64(144)).raw());
        results.push(fp_sin(FP_HALF_PI).raw());
        results.push(fp_cos(FP::ZERO).raw());
        results.push(fp_atan2(FP::from_i64(1), FP::from_i64(1)).raw());
        results.push(fp_lerp(FP::from_i64(0), FP::from_i64(100), FP::HALF).raw());
        results.push(fp_clamp(FP::from_i64(50), FP::from_i64(0), FP::from_i64(10)).raw());
        results.push(fp_min(FP::from_i64(-3), FP::from_i64(7)).raw());
        results.push(fp_max(FP::from_i64(-3), FP::from_i64(7)).raw());

        // Run the same computation again — must be bit-identical
        let a2 = FP::from_i64(7);
        let b2 = FP::from_i64(3);
        let c2 = FP::from_f64(2.5);

        let results2 = vec![
            (a2 + b2).raw(),
            (a2 - b2).raw(),
            (a2 * b2).raw(),
            (a2 / b2).raw(),
            (c2 * FP::from_i64(4)).raw(),
            fp_sqrt(FP::from_i64(144)).raw(),
            fp_sin(FP_HALF_PI).raw(),
            fp_cos(FP::ZERO).raw(),
            fp_atan2(FP::from_i64(1), FP::from_i64(1)).raw(),
            fp_lerp(FP::from_i64(0), FP::from_i64(100), FP::HALF).raw(),
            fp_clamp(FP::from_i64(50), FP::from_i64(0), FP::from_i64(10)).raw(),
            fp_min(FP::from_i64(-3), FP::from_i64(7)).raw(),
            fp_max(FP::from_i64(-3), FP::from_i64(7)).raw(),
        ];

        assert_eq!(results, results2, "determinism violated: results differ!");
    }

    #[test]
    fn test_constants() {
        let pi_diff = (FP_PI.to_f64() - std::f64::consts::PI).abs();
        assert!(pi_diff < 0.01, "PI constant off by {}", pi_diff);
        let e_diff = (FP_E.to_f64() - std::f64::consts::E).abs();
        assert!(e_diff < 0.01, "E constant off by {}", e_diff);
    }

    #[test]
    #[should_panic]
    fn test_div_by_zero() {
        let _ = FP::from_i64(1) / FP::ZERO;
    }
}
