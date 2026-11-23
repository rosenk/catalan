package quadratic

import (
	"errors"
	"math/big"

	"github.com/rosenk/catalan/rings"
)

// Грешки от AutoSolve.
var (
	ErrNoRealRoots = errors.New("уравнението няма реални корени (D < 0)")
	ErrDegenerate  = errors.New("c1 и c2 са 0, уравнението няма неизвестно")
	ErrBadDigits   = errors.New("броят знаци трябва да е поне 1")
	errNoProgress  = errors.New("вътрешна грешка: итерацията не е сходяща")
)

const (
	// StepTerms е броят членове на реда на Каталан в една стъпка.
	// Стъпка с K члена е метод от ред K+1.
	StepTerms = 8
	maxSteps  = 64 // предпазна граница, на практика стигат 3–5 стъпки
	guardBits = 32 // допълнителни битове над исканата точност
)

// Discriminant връща D = c1² - 4*c0*c2 = c1²(1 - 4t). D не се променя при Shift.
func Discriminant(c0, c1, c2 *big.Rat) *big.Rat {
	r := rings.BigRat{}
	return r.Sub(r.Mul(c1, c1), r.Mul(big.NewRat(4, 1), r.Mul(c0, c2)))
}

// Shift замества x = y + s в c0 - c1*x + c2*x² = 0
// и връща коефициентите на уравнението за y в същата форма.
func Shift(c0, c1, c2, s *big.Rat) (d0, d1, d2 *big.Rat) {
	r := rings.BigRat{}
	d0 = r.Add(r.Sub(c0, r.Mul(c1, s)), r.Mul(c2, r.Mul(s, s)))
	d1 = r.Sub(c1, r.Mul(big.NewRat(2, 1), r.Mul(c2, s)))
	d2 = new(big.Rat).Set(c2)
	return d0, d1, d2
}

// Result описва решението от AutoSolve.
type Result struct {
	Start      *big.Rat // началната точка от StartShift
	StartParam *big.Rat // параметърът на реда t в началната точка, винаги в [0, 5/36)
	Steps      int      // брой стъпки с StepTerms члена
	X1         *big.Rat // коренът от реда на Каталан (по-далечният от 0)
	X2         *big.Rat // вторият корен по Виет, nil при линейно уравнение
	Double     bool     // D = 0, двоен корен
	Linear     bool     // c2 = 0
}

// roundSig закръгля x до около bits значещи бита със знаменател степен на 2.
// Така числителят и знаменателят не растат повече от нужното.
func roundSig(x *big.Rat, bits int) *big.Rat {
	if x.Sign() == 0 {
		return new(big.Rat)
	}
	shift := bits - (x.Num().BitLen() - x.Denom().BitLen())
	num := new(big.Int).Abs(x.Num())
	den := new(big.Int).Set(x.Denom())
	if shift >= 0 {
		num.Lsh(num, uint(shift))
	} else {
		den.Lsh(den, uint(-shift))
	}
	q, m := new(big.Int).QuoRem(num, den, new(big.Int))
	if m.Lsh(m, 1).Cmp(den) >= 0 {
		q.Add(q, big.NewInt(1))
	}
	if x.Sign() < 0 {
		q.Neg(q)
	}
	res := new(big.Rat).SetInt(q)
	if shift >= 0 {
		return res.Quo(res, new(big.Rat).SetInt(new(big.Int).Lsh(big.NewInt(1), uint(shift))))
	}
	return res.Mul(res, new(big.Rat).SetInt(new(big.Int).Lsh(big.NewInt(1), uint(-shift))))
}

// floorDiv2 връща ⌊e/2⌋ и за отрицателни e.
func floorDiv2(e int) int {
	if e < 0 {
		return -((-e + 1) / 2)
	}
	return e / 2
}

func abs(n int) int {
	if n < 0 {
		return -n
	}
	return n
}

// StartShift връща начална точка s = v ± h, където v = c1/(2*c2) е средата между корените.
// Знакът е знакът на v (плюс при v = 0), затова s е до корена, който е по-далеч от 0.
//
// Нека H е разстоянието от v до корените. Тогава R = H² = D/(4*c2²) е рационално число.
// Програмата не пресмята H. Тя взима k = (битове на числителя) - (битове на знаменателя) на R
// и m = ⌊k/2⌋. Тогава 2^(k-1) < R < 2^(k+1).
//   - При четно k: H е между 2^m/√2 и 2^m·√2. Кандидатите за h са 2^m и (3/2)·2^m.
//   - При нечетно k: H е между 2^m и 2^(m+1). Кандидатите за h са (3/2)·2^m и 2^(m+1).
//
// Програмата взима първия кандидат, при който s е извън корените. Това е тест на знака:
// t(s) = p(s)*c2 / p'(s)² >= 0. И в двата случая H <= h < (3/2)·H.
// Нека ρ = (h - H)/(h + H) е отношението на разстоянията от s до двата корена. Тогава ρ < 1/5,
// и параметърът на реда t = ρ/(1+ρ)² е в [0, 5/36).
//
// При t >= 0 стъпките в AutoSolve пазят t >= 0, защото p(s + y) = p(s)·E_M(t),
// а E_M(t) = 1 - P_M(t) + t·P_M(t)² има само положителни коефициенти.
// Изисква c2 != 0 и D > 0.
func StartShift(c0, c1, c2 *big.Rat) *big.Rat {
	D := Discriminant(c0, c1, c2)
	R := new(big.Rat).Quo(D, new(big.Rat).Mul(big.NewRat(4, 1), new(big.Rat).Mul(c2, c2)))
	k := R.Num().BitLen() - R.Denom().BitLen()
	m := floorDiv2(k)

	pow := new(big.Rat).SetInt(new(big.Int).Lsh(big.NewInt(1), uint(abs(m))))
	if m < 0 {
		pow.Inv(pow)
	}
	threeHalves := new(big.Rat).Mul(big.NewRat(3, 2), pow)
	first, second := pow, threeHalves
	if k%2 != 0 {
		first, second = threeHalves, new(big.Rat).Mul(big.NewRat(2, 1), pow)
	}

	v := new(big.Rat).Quo(c1, new(big.Rat).Mul(big.NewRat(2, 1), c2))
	at := func(h *big.Rat) *big.Rat {
		if v.Sign() < 0 {
			return new(big.Rat).Sub(v, h)
		}
		return new(big.Rat).Add(v, h)
	}
	s := at(first)
	if d0, _, _ := Shift(c0, c1, c2, s); new(big.Rat).Mul(d0, c2).Sign() < 0 {
		s = at(second)
	}
	return s
}

// AutoSolve решава c0 - c1*x + c2*x² = 0 с относителна точност около digits десетични знака.
// Програмата използва само реда на Каталан и аритметика в ℚ, без корен квадратен:
//  1. Ако D < 0, корени няма. Ако D = 0, α(1/4) = 2 и x = 2*c0/c1 = c1/(2*c2).
//  2. Програмата започва от StartShift. Там 0 <= t < 5/36 и t остава >= 0 на всяка стъпка.
//  3. На всяка стъпка тя измества началото в s, сумира StepTerms члена от реда за y
//     и премества s в s + y. Стъпката е метод от ред StepTerms+1.
//     От t = 5/36 първата стъпка дава около 3 знака, втората около 30, третата около 270.
//  4. Програмата спира, когато |y| е под |s|·2^-bits.
//  5. Вторият корен идва от формулата на Виет x1*x2 = c0/c2.
//     Тя пази относителната точност и на малкия по модул корен.
func AutoSolve(c0, c1, c2 *big.Rat, digits int) (*Result, error) {
	if digits < 1 {
		return nil, ErrBadDigits
	}
	if c2.Sign() == 0 {
		if c1.Sign() == 0 {
			return nil, ErrDegenerate
		}
		return &Result{X1: new(big.Rat).Quo(c0, c1), Linear: true}, nil
	}

	D := Discriminant(c0, c1, c2)
	switch D.Sign() {
	case -1:
		return nil, ErrNoRealRoots
	case 0:
		x := new(big.Rat).Quo(c1, new(big.Rat).Mul(big.NewRat(2, 1), c2))
		return &Result{X1: x, X2: x, Double: true}, nil
	}

	bits := digits*10/3 + 1 + guardBits // log2(10) < 10/3

	s := StartShift(c0, c1, c2)
	res := &Result{Start: new(big.Rat).Set(s)}
	d0, d1, d2 := Shift(c0, c1, c2, s)
	res.StartParam = Param(T(d0, d1, d2))

	// prec = |s| · 2^-bits
	prec := new(big.Rat).Abs(s)
	prec.Quo(prec, new(big.Rat).SetInt(new(big.Int).Lsh(big.NewInt(1), uint(bits))))

	for {
		if !Converges(d0, d1, d2) {
			return nil, errNoProgress
		}
		y := roundSig(Solve(d0, d1, d2, StepTerms), bits)
		s.Add(s, y)
		res.Steps++
		if new(big.Rat).Abs(y).Cmp(prec) <= 0 {
			break
		}
		if res.Steps >= maxSteps {
			return nil, errNoProgress
		}
		d0, d1, d2 = Shift(c0, c1, c2, s)
	}

	res.X1 = s
	res.X2 = new(big.Rat).Quo(c0, new(big.Rat).Mul(c2, s))
	return res, nil
}
