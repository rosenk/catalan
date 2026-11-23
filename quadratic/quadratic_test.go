package quadratic_test

import (
	"math/big"
	"math/rand"
	"testing"

	"github.com/rosenk/catalan/quadratic"
)

func TestSolveConverges(t *testing.T) {
	c0, c1, c2 := big.NewRat(1, 1), big.NewRat(5, 1), big.NewRat(1, 1)
	prev := new(big.Rat).SetInt64(1)
	for _, terms := range []int{5, 10, 20, 40} {
		x := quadratic.Solve(c0, c1, c2, terms)
		res := quadratic.Residual(c0, c1, c2, x)
		if res.Cmp(prev) >= 0 {
			t.Fatalf("terms=%d: residual %s did not decrease", terms, res.FloatString(25))
		}
		prev = res
	}
	if prev.Cmp(big.NewRat(1, 1_000_000_000_000_000)) >= 0 {
		t.Errorf("residual after 40 terms too large: %s", prev.FloatString(25))
	}
}

func TestConverges(t *testing.T) {
	// 1 - 5x + x² = 0: t = 1/25, сходящ
	if !quadratic.Converges(big.NewRat(1, 1), big.NewRat(5, 1), big.NewRat(1, 1)) {
		t.Error("1 - 5x + x²: expected convergence")
	}
	// x² + 5x - 7 = 0 → c0=-7, c1=-5, c2=1: t = -7/25, сходящ чрез u = 7/53
	if !quadratic.Converges(big.NewRat(-7, 1), big.NewRat(-5, 1), big.NewRat(1, 1)) {
		t.Error("x² + 5x - 7: expected convergence")
	}
	// x² + 1 = 0 с изместване 1 → t > 1/4, разходящ
	if quadratic.Converges(big.NewRat(2, 1), big.NewRat(2, 1), big.NewRat(1, 1)) {
		t.Error("t = 1/2: expected divergence")
	}
	// c1 = 0
	if quadratic.Converges(big.NewRat(1, 1), new(big.Rat), big.NewRat(1, 1)) {
		t.Error("c1 = 0: expected false")
	}
}

// relRes връща |p(x)| / (|c0| + |c1*x| + |c2*x²|).
func relRes(c0, c1, c2, x *big.Rat) *big.Rat {
	num := quadratic.Residual(c0, c1, c2, x)
	den := new(big.Rat).Abs(c0)
	den.Add(den, new(big.Rat).Abs(new(big.Rat).Mul(c1, x)))
	den.Add(den, new(big.Rat).Abs(new(big.Rat).Mul(c2, new(big.Rat).Mul(x, x))))
	return num.Quo(num, den)
}

func pow10(n int) *big.Rat {
	p := new(big.Rat).SetInt(new(big.Int).Exp(big.NewInt(10), big.NewInt(int64(abs(n))), nil))
	if n < 0 {
		p.Inv(p)
	}
	return p
}

func TestAutoSolve(t *testing.T) {
	tol := pow10(-30)
	cases := []struct {
		name       string
		c0, c1, c2 *big.Rat
	}{
		{"x²+5x-7", big.NewRat(-7, 1), big.NewRat(-5, 1), big.NewRat(1, 1)},
		{"1-5x+x²", big.NewRat(1, 1), big.NewRat(5, 1), big.NewRat(1, 1)},
		{"x²-2 (c1=0)", big.NewRat(-2, 1), new(big.Rat), big.NewRat(1, 1)},
		{"близки корени", big.NewRat(9999, 1), big.NewRat(200, 1), big.NewRat(1, 1)},
		{"далечни корени", big.NewRat(-1_000_000, 1), big.NewRat(1, 1), big.NewRat(1, 1)},
		{"дробни", big.NewRat(-1, 3), big.NewRat(7, 2), big.NewRat(5, 7)},
		{"t = -100", big.NewRat(-100, 1), big.NewRat(1, 1), big.NewRat(1, 1)},
		{"x²-10^20x+1", big.NewRat(1, 1), new(big.Rat).SetInt(new(big.Int).Exp(big.NewInt(10), big.NewInt(20), nil)), big.NewRat(1, 1)},
		{"c0 = 0", new(big.Rat), big.NewRat(3, 1), big.NewRat(1, 1)},
		{"корени около ±2^250", big.NewRat(-1, 1), new(big.Rat), new(big.Rat).SetFrac(big.NewInt(1), new(big.Int).Lsh(big.NewInt(1), 500))},
		{"корени около ±2^-250", new(big.Rat).SetFrac(big.NewInt(-1), new(big.Int).Lsh(big.NewInt(1), 500)), new(big.Rat), big.NewRat(1, 1)},
	}
	for _, c := range cases {
		res, err := quadratic.AutoSolve(c.c0, c.c1, c.c2, 30)
		if err != nil {
			t.Errorf("%s: %v", c.name, err)
			continue
		}
		for _, x := range []*big.Rat{res.X1, res.X2} {
			if x.Sign() == 0 && c.c0.Sign() == 0 {
				continue
			}
			if r := relRes(c.c0, c.c1, c.c2, x); r.Cmp(tol) >= 0 {
				t.Errorf("%s: x=%s relative residual %s", c.name, x.FloatString(10), r.FloatString(40))
			}
		}
	}
}

func TestAutoSolveErrors(t *testing.T) {
	// x² + 1 = 0
	if _, err := quadratic.AutoSolve(big.NewRat(1, 1), new(big.Rat), big.NewRat(1, 1), 20); err != quadratic.ErrNoRealRoots {
		t.Errorf("x²+1: err = %v", err)
	}
	// x² - 2x + 1 = 0 → двоен корен 1
	res, err := quadratic.AutoSolve(big.NewRat(1, 1), big.NewRat(2, 1), big.NewRat(1, 1), 20)
	if err != nil || !res.Double || res.X1.Cmp(big.NewRat(1, 1)) != 0 {
		t.Errorf("(x-1)²: res=%+v err=%v", res, err)
	}
	// 1 - x = 0 → линейно, x = 1
	res, err = quadratic.AutoSolve(big.NewRat(1, 1), big.NewRat(1, 1), new(big.Rat), 20)
	if err != nil || !res.Linear || res.X1.Cmp(big.NewRat(1, 1)) != 0 {
		t.Errorf("линейно: res=%+v err=%v", res, err)
	}
	if _, err := quadratic.AutoSolve(big.NewRat(-2, 1), new(big.Rat), big.NewRat(1, 1), 0); err != quadratic.ErrBadDigits {
		t.Errorf("digits=0: err = %v", err)
	}
}

func TestMobius(t *testing.T) {
	// x² + 5x - 7 = 0 без изместване: u = 7/53, грешката е около 0.528^n
	c0, c1, c2 := big.NewRat(-7, 1), big.NewRat(-5, 1), big.NewRat(1, 1)
	if u := quadratic.Param(quadratic.T(c0, c1, c2)); u.Cmp(big.NewRat(7, 53)) != 0 {
		t.Fatalf("u = %s, want 7/53", u.RatString())
	}
	x := quadratic.Solve(c0, c1, c2, 80)
	if r := quadratic.Residual(c0, c1, c2, x); r.Cmp(big.NewRat(1, 1_000_000_000_000_000)) >= 0 {
		t.Errorf("x = %s, residual %s", x.FloatString(20), r.FloatString(25))
	}
}

// randRat връща случайна ненулева дроб от порядък 10^exp.
func randRat(rng *rand.Rand, exp int) *big.Rat {
	num := big.NewInt(rng.Int63n(1_000_000) + 1)
	if rng.Intn(2) == 0 {
		num.Neg(num)
	}
	r := new(big.Rat).SetFrac(num, big.NewInt(rng.Int63n(1_000_000)+1))
	return r.Mul(r, pow10(exp))
}

func abs(n int) int {
	if n < 0 {
		return -n
	}
	return n
}

// TestStartGuarantee решава случайни уравнения a(x - r1)(x - r2) = 0 с корени от 10^-50 до 10^50
// и разстояние между тях до 10^-300. Тестът проверява границата [0, 5/36), броя стъпки и точността.
func TestStartGuarantee(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	bound := big.NewRat(5, 36)
	tol := pow10(-30)
	for i := 0; i < 300; i++ {
		a := randRat(rng, 0)
		r1 := randRat(rng, rng.Intn(101)-50)
		r2 := new(big.Rat).Add(r1, randRat(rng, rng.Intn(351)-300))
		c2 := a
		c1 := new(big.Rat).Mul(a, new(big.Rat).Add(r1, r2))
		c0 := new(big.Rat).Mul(a, new(big.Rat).Mul(r1, r2))

		res, err := quadratic.AutoSolve(c0, c1, c2, 30)
		if err != nil {
			t.Fatalf("case %d: %v", i, err)
		}
		if res.StartParam.Sign() < 0 || res.StartParam.Cmp(bound) >= 0 {
			t.Fatalf("case %d: start t = %s, want [0, 5/36)", i, res.StartParam.FloatString(6))
		}
		if res.Steps > 5 {
			t.Fatalf("case %d: %d steps", i, res.Steps)
		}
		for _, x := range []*big.Rat{res.X1, res.X2} {
			if r := relRes(c0, c1, c2, x); r.Cmp(tol) >= 0 {
				t.Fatalf("case %d: relative residual %s", i, r.FloatString(40))
			}
		}
	}
}

func TestVeryCloseRoots(t *testing.T) {
	// (x-1)² = 10^-2000, корени 1 ± 10^-1000
	c0 := new(big.Rat).Sub(big.NewRat(1, 1), pow10(-2000))
	res, err := quadratic.AutoSolve(c0, big.NewRat(2, 1), big.NewRat(1, 1), 1030)
	if err != nil {
		t.Fatal(err)
	}
	if res.Steps > 5 {
		t.Errorf("steps = %d", res.Steps)
	}
	// x1 - 1 и 1 - x2 трябва да са 10^-1000 с относителна грешка под 10^-20
	for _, d := range []*big.Rat{
		new(big.Rat).Sub(res.X1, big.NewRat(1, 1)),
		new(big.Rat).Sub(big.NewRat(1, 1), res.X2),
	} {
		rel := new(big.Rat).Quo(new(big.Rat).Sub(d, pow10(-1000)), pow10(-1000))
		if rel.Abs(rel).Cmp(pow10(-20)) >= 0 {
			t.Errorf("relative error %s", rel.FloatString(30))
		}
	}
}
