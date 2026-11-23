package fps_test

import (
	"math/big"
	"testing"

	"github.com/rosenk/catalan/fps"
	"github.com/rosenk/catalan/rings"
)

func TestCatalanNumbers(t *testing.T) {
	want := []int64{1, 1, 2, 5, 14, 42, 132, 429, 1430, 4862, 16796}
	alpha := fps.Catalan[*big.Int](rings.BigInt{})
	for n, w := range want {
		if got := alpha.Coeff(n); got.Cmp(big.NewInt(w)) != 0 {
			t.Errorf("C_%d = %s, want %d", n, got, w)
		}
	}
}

func TestIdentityZ(t *testing.T) {
	r := rings.BigInt{}
	alpha := fps.Catalan[*big.Int](r)
	eq := fps.One[*big.Int](r).Sub(alpha).Add(fps.T[*big.Int](r).Mul(alpha.Mul(alpha)))
	for n := 0; n < 30; n++ {
		if c := eq.Coeff(n); c.Sign() != 0 {
			t.Fatalf("coefficient %d = %s, want 0", n, c)
		}
	}
}

func TestIdentityGF2(t *testing.T) {
	r := rings.GF2Ring{}
	alpha := fps.Catalan[rings.GF2](r)
	eq := fps.One[rings.GF2](r).Sub(alpha).Add(fps.T[rings.GF2](r).Mul(alpha.Mul(alpha)))
	for n := 0; n < 30; n++ {
		if c := eq.Coeff(n); c != 0 {
			t.Fatalf("coefficient %d = %d, want 0", n, c)
		}
	}
}

func TestEvalHorner(t *testing.T) {
	r := rings.BigInt{}
	// (1 + t)^2 = 1 + 2t + t^2; при t = 3 → 16
	p := fps.One[*big.Int](r).Add(fps.T[*big.Int](r))
	if got := p.Mul(p).Eval(big.NewInt(3), 5); got.Cmp(big.NewInt(16)) != 0 {
		t.Errorf("Eval = %s, want 16", got)
	}
}
