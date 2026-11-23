// Package quadratic решава c0 - c1*x + c2*x^2 = 0 изцяло в ℚ, без корен квадратен.
package quadratic

import (
	"math/big"

	"github.com/rosenk/catalan/fps"
	"github.com/rosenk/catalan/rings"
)

var quarter = big.NewRat(1, 4)

// T връща t = c0*c2 / c1². c1 не трябва да е 0.
func T(c0, c1, c2 *big.Rat) *big.Rat {
	return new(big.Rat).Quo(new(big.Rat).Mul(c0, c2), new(big.Rat).Mul(c1, c1))
}

// Param връща числото, в което Solve сумира реда на Каталан:
// t при t >= 0 и u = -t / (1 - 4t) при t < 0 (замяна на Мьобиус).
// За всяко t < 0 важи 0 < u < 1/4.
func Param(t *big.Rat) *big.Rat {
	if t.Sign() >= 0 {
		return new(big.Rat).Set(t)
	}
	oneMinus4t := new(big.Rat).Sub(big.NewRat(1, 1), new(big.Rat).Mul(big.NewRat(4, 1), t))
	return new(big.Rat).Quo(new(big.Rat).Neg(t), oneMinus4t)
}

// Converges казва дали Solve е сходяща: c1 != 0 и t < 1/4.
// Това е същото като c1 != 0 и D > 0.
func Converges(c0, c1, c2 *big.Rat) bool {
	return c1.Sign() != 0 && T(c0, c1, c2).Cmp(quarter) < 0
}

// Solve връща приближение на корена на c0 - c1*x + c2*x² = 0, който е най-близо до 0,
// чрез първите terms члена на реда на Каталан: x = (c0/c1) * α(t).
//
// При 0 <= t < 1/4 редът α(t) е сходящ директно.
// При t < 0 Solve сумира β = α(u) за u = -t/(1-4t) и връща α = β / (1 - 4t + 2tβ).
// Грешката намалява като (4*Param(t))^terms.
// Извикващият трябва да е проверил Converges.
func Solve(c0, c1, c2 *big.Rat, terms int) *big.Rat {
	r := rings.BigRat{}
	t := T(c0, c1, c2)
	cat := fps.Catalan[*big.Rat](r)

	var alpha *big.Rat
	if t.Sign() >= 0 {
		alpha = cat.Eval(t, terms)
	} else {
		oneMinus4t := r.Sub(r.One(), r.Mul(big.NewRat(4, 1), t))
		beta := cat.Eval(Param(t), terms)
		denom := r.Add(oneMinus4t, r.Mul(big.NewRat(2, 1), r.Mul(t, beta)))
		alpha = r.Mul(beta, r.Add(r.Zero(), new(big.Rat).Inv(denom)))
	}
	return r.Mul(new(big.Rat).Quo(c0, c1), alpha)
}

// Residual връща |c0 - c1*x + c2*x²|.
func Residual(c0, c1, c2, x *big.Rat) *big.Rat {
	r := rings.BigRat{}
	res := r.Add(r.Sub(c0, r.Mul(c1, x)), r.Mul(c2, r.Mul(x, x)))
	return res.Abs(res)
}
