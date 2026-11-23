package fps

import "math/big"

// Catalan генерира потока C_n по рекурентността C_{m+1} = Σ_{k=0..m} C_k * C_{m-k}.
// Числата се пресмятат в ℤ и се кешират, после се изобразяват в пръстена R.
// Върнатият ред не е безопасен за едновременна употреба от няколко goroutine.
func Catalan[R any](r Ring[R]) FPS[R] {
	cache := []*big.Int{big.NewInt(1)}

	catalanBig := func(n int) *big.Int {
		for len(cache) <= n {
			m := len(cache) - 1
			sum := new(big.Int)
			prod := new(big.Int)
			for k := 0; k <= m; k++ {
				prod.Mul(cache[k], cache[m-k])
				sum.Add(sum, prod)
			}
			cache = append(cache, sum)
		}
		return cache[n]
	}

	return New(r, func(n int) R { return r.FromBigInt(catalanBig(n)) })
}
