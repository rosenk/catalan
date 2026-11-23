// Package fps реализира формални степенни редове над произволен комутативен пръстен.
package fps

import "math/big"

// Ring е минимален комутативен пръстен (еквивалент на `class CommRing R`).
type Ring[R any] interface {
	Zero() R
	One() R
	Add(a, b R) R
	Sub(a, b R) R
	Mul(a, b R) R
	FromBigInt(n *big.Int) R
}

// FPS е формален степенен ред (еквивалент на `structure FPS R`).
type FPS[R any] struct {
	ring  Ring[R]
	coeff func(n int) R
}

// New създава ред от пръстен и функция за коефициентите.
func New[R any](r Ring[R], coeff func(n int) R) FPS[R] {
	return FPS[R]{ring: r, coeff: coeff}
}

// Coeff връща n-тия коефициент.
func (f FPS[R]) Coeff(n int) R {
	return f.coeff(n)
}

// Add връща f + g.
func (f FPS[R]) Add(g FPS[R]) FPS[R] {
	r := f.ring
	return New(r, func(n int) R { return r.Add(f.coeff(n), g.coeff(n)) })
}

// Sub връща f - g.
func (f FPS[R]) Sub(g FPS[R]) FPS[R] {
	r := f.ring
	return New(r, func(n int) R { return r.Sub(f.coeff(n), g.coeff(n)) })
}

// Mul е произведението на Коши (еквивалент на `sumFin (n + 1) ...`).
func (f FPS[R]) Mul(g FPS[R]) FPS[R] {
	r := f.ring
	return New(r, func(n int) R {
		sum := r.Zero()
		for k := 0; k <= n; k++ {
			sum = r.Add(sum, r.Mul(f.coeff(k), g.coeff(n-k)))
		}
		return sum
	})
}

// Eval пресмята частичната сума Σ_{k=0}^{terms-1} a_k * t^k по схемата на Хорнер.
func (f FPS[R]) Eval(tVal R, terms int) R {
	r := f.ring
	acc := r.Zero()
	for k := terms - 1; k >= 0; k-- {
		acc = r.Add(r.Mul(acc, tVal), f.coeff(k))
	}
	return acc
}

// One е константният ред 1.
func One[R any](r Ring[R]) FPS[R] {
	return New(r, func(n int) R {
		if n == 0 {
			return r.One()
		}
		return r.Zero()
	})
}

// T е редът на променливата t.
func T[R any](r Ring[R]) FPS[R] {
	return New(r, func(n int) R {
		if n == 1 {
			return r.One()
		}
		return r.Zero()
	})
}
