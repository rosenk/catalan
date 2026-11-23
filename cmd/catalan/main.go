// Команда catalan демонстрира реда на Каталан над ℤ, 𝔽₂ и ℚ
// и решава уравнения от вида a*x² + b*x + c = 0.
package main

import (
	"flag"
	"fmt"
	"math/big"
	"os"

	"github.com/rosenk/catalan/fps"
	"github.com/rosenk/catalan/quadratic"
	"github.com/rosenk/catalan/rings"
)

// identity връща 1 - α + t·α² за реда на Каталан α над пръстена r.
func identity[R any](r fps.Ring[R]) (alpha, eq fps.FPS[R]) {
	alpha = fps.Catalan(r)
	eq = fps.One(r).Sub(alpha).Add(fps.T(r).Mul(alpha.Mul(alpha)))
	return alpha, eq
}

func parseRat(name, s string) *big.Rat {
	v, ok := new(big.Rat).SetString(s)
	if !ok {
		fmt.Fprintf(os.Stderr, "невалидно число за -%s: %q\n", name, s)
		os.Exit(2)
	}
	return v
}

func demo(checkN int) {
	alphaZ, eqZ := identity[*big.Int](rings.BigInt{})

	fmt.Print("Catalan C_0..C_9: ")
	for n := 0; n < 10; n++ {
		fmt.Printf("%s ", alphaZ.Coeff(n))
	}
	fmt.Println()

	allZeroZ := true
	for n := 0; n < checkN; n++ {
		if eqZ.Coeff(n).Sign() != 0 {
			allZeroZ = false
		}
	}
	fmt.Printf("1 - α + t·α² == 0 в ℤ[[t]] за n=0..%d: %v\n", checkN-1, allZeroZ)

	_, eq2 := identity[rings.GF2](rings.GF2Ring{})
	allZeroGF2 := true
	for n := 0; n < checkN; n++ {
		if eq2.Coeff(n) != 0 {
			allZeroGF2 = false
		}
	}
	fmt.Printf("1 - α + t·α² == 0 в 𝔽₂[[t]] за n=0..%d: %v\n", checkN-1, allZeroGF2)
}

// format показва x с digits значещи цифри.
func format(x *big.Rat, digits int) string {
	if x.Sign() == 0 {
		return "0"
	}
	f := new(big.Float).SetPrec(uint(digits*4 + 64)).SetRat(x)
	return f.Text('g', digits)
}

// relResidual връща |c0 - c1*x + c2*x²| / (|c0| + |c1*x| + |c2*x²|).
func relResidual(c0, c1, c2, x *big.Rat) string {
	num := quadratic.Residual(c0, c1, c2, x)
	if num.Sign() == 0 {
		return "0"
	}
	den := new(big.Rat).Abs(c0)
	den.Add(den, new(big.Rat).Abs(new(big.Rat).Mul(c1, x)))
	den.Add(den, new(big.Rat).Abs(new(big.Rat).Mul(c2, new(big.Rat).Mul(x, x))))
	return new(big.Float).SetPrec(64).SetRat(num.Quo(num, den)).Text('e', 2)
}

func main() {
	checkN := flag.Int("n", 20, "брой коефициенти за проверка на тъждеството")
	digits := flag.Int("digits", 30, "брой значещи цифри на корените")
	as := flag.String("a", "1", "коефициент a в a*x² + b*x + c = 0 (цяло, дроб 7/3 или десетично)")
	bs := flag.String("b", "-5", "коефициент b в a*x² + b*x + c = 0")
	cs := flag.String("c", "1", "коефициент c в a*x² + b*x + c = 0")
	flag.Parse()

	if *digits < 1 {
		fmt.Fprintln(os.Stderr, "-digits трябва да е поне 1")
		os.Exit(2)
	}

	userEq := false
	flag.Visit(func(f *flag.Flag) {
		switch f.Name {
		case "a", "b", "c":
			userEq = true
		}
	})
	if !userEq {
		demo(*checkN)
	}

	a, b, c := parseRat("a", *as), parseRat("b", *bs), parseRat("c", *cs)
	fmt.Printf("уравнение: (%s)x² + (%s)x + (%s) = 0\n", a.RatString(), b.RatString(), c.RatString())

	// Вътрешната форма е c0 - c1*x + c2*x² = 0, тоест c0 = c, c1 = -b, c2 = a.
	c0, c1, c2 := c, new(big.Rat).Neg(b), a

	res, err := quadratic.AutoSolve(c0, c1, c2, *digits)
	if err != nil {
		fmt.Fprintln(os.Stderr, "методът не може да реши уравнението:", err)
		os.Exit(1)
	}

	switch {
	case res.Linear:
		fmt.Printf("линейно уравнение, x = %s\n", res.X1.RatString())
		return
	case res.Double:
		fmt.Printf("D = 0, двоен корен x = %s\n", res.X1.RatString())
		return
	}

	sp, _ := res.StartParam.Float64()
	fmt.Printf("начална точка s = %s, параметър на реда там ≈ %.4g\n", format(res.Start, 12), sp)
	fmt.Printf("стъпки по %d члена от реда на Каталан: %d\n", quadratic.StepTerms, res.Steps)
	fmt.Printf("x1 = %s\n", format(res.X1, *digits))
	fmt.Printf("   относителен остатък: %s\n", relResidual(c0, c1, c2, res.X1))
	fmt.Printf("x2 = %s (по Виет, x1·x2 = c/a)\n", format(res.X2, *digits))
	fmt.Printf("   относителен остатък: %s\n", relResidual(c0, c1, c2, res.X2))
}
