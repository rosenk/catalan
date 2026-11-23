// Package rings съдържа конкретни инстанции на fps.Ring: ℤ, ℚ и 𝔽₂.
package rings

import (
	"math/big"

	"github.com/rosenk/catalan/fps"
)

// Проверки по време на компилация, че типовете изпълняват интерфейса.
var (
	_ fps.Ring[*big.Int] = BigInt{}
	_ fps.Ring[*big.Rat] = BigRat{}
	_ fps.Ring[GF2]      = GF2Ring{}
)

// BigInt е пръстенът ℤ.
type BigInt struct{}

func (BigInt) Zero() *big.Int                 { return big.NewInt(0) }
func (BigInt) One() *big.Int                  { return big.NewInt(1) }
func (BigInt) Add(a, b *big.Int) *big.Int     { return new(big.Int).Add(a, b) }
func (BigInt) Sub(a, b *big.Int) *big.Int     { return new(big.Int).Sub(a, b) }
func (BigInt) Mul(a, b *big.Int) *big.Int     { return new(big.Int).Mul(a, b) }
func (BigInt) FromBigInt(n *big.Int) *big.Int { return new(big.Int).Set(n) }

// BigRat е полето ℚ.
type BigRat struct{}

func (BigRat) Zero() *big.Rat                 { return new(big.Rat) }
func (BigRat) One() *big.Rat                  { return big.NewRat(1, 1) }
func (BigRat) Add(a, b *big.Rat) *big.Rat     { return new(big.Rat).Add(a, b) }
func (BigRat) Sub(a, b *big.Rat) *big.Rat     { return new(big.Rat).Sub(a, b) }
func (BigRat) Mul(a, b *big.Rat) *big.Rat     { return new(big.Rat).Mul(a, b) }
func (BigRat) FromBigInt(n *big.Int) *big.Rat { return new(big.Rat).SetInt(n) }

// GF2 е елемент на крайното поле с характеристика 2
// (където класическата формула дели на 0).
type GF2 uint8

// GF2Ring е полето 𝔽₂.
type GF2Ring struct{}

func (GF2Ring) Zero() GF2                 { return 0 }
func (GF2Ring) One() GF2                  { return 1 }
func (GF2Ring) Add(a, b GF2) GF2          { return (a ^ b) & 1 }
func (GF2Ring) Sub(a, b GF2) GF2          { return (a ^ b) & 1 }
func (GF2Ring) Mul(a, b GF2) GF2          { return (a & b) & 1 }
func (GF2Ring) FromBigInt(n *big.Int) GF2 { return GF2(n.Bit(0)) }
