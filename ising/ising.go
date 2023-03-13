package ising

import (
	"math"
	"math/rand"
)

type ising struct {
	size   int
	t      float64
	v      []int
	lookup []float64
	nn     []int
}

func New(size int, t float64) *ising {
	l := ising{size: size, t: t}
	l.v = make([]int, size*size)
	for i := range l.v {
		l.v[i] = 1
	}

	l.nn = make([]int, 4*size*size)
	calc_nn(l.nn, size)
	l.lookup = make([]float64, 2)
	calc_lookup(l.lookup, t)

	return &l
}

func (l *ising) Set_t(t float64) {
	l.t = t
	calc_lookup(l.lookup, t)
}
func (l ising) Get_t() float64 {
	return l.t
}
func (l ising) Get_Size() int {
	return l.size
}

func calc_nn(nn []int, size int) {
	for i := 0; i < size*size; i++ {
		//left
		if i%size == 0 {
			nn[4*i] = i + (size - 1)
		} else {
			nn[4*i] = i - 1
		}
		// bottom
		if i >= size*(size-1) {
			nn[4*i+1] = i % size
		} else {
			nn[4*i+1] = i + size
		}
		// right
		if i%size == size-1 {
			nn[4*i+2] = i - (size - 1)
		} else {
			nn[4*i+2] = i + 1
		}
		//top
		if i < size {
			nn[4*i+3] = size*(size-1) + i
		} else {
			nn[4*i+3] = i - size
		}
	}
}

func calc_lookup(lookup []float64, t float64) {
	lookup[0] = math.Exp(-4.0 / t)
	lookup[1] = math.Exp(-8.0 / t)
}

func (l *ising) dE(i int) int {
	res := 2 * l.v[i] * l.v[l.nn[4*i]]
	res += 2 * l.v[i] * l.v[l.nn[4*i+1]]
	res += 2 * l.v[i] * l.v[l.nn[4*i+2]]
	res += 2 * l.v[i] * l.v[l.nn[4*i+3]]
	return res
}

func (l *ising) Sweep(r *rand.Rand) {
	for i := 0; i < l.size*l.size; i++ {
		de := l.dE(i)
		if de <= 0 {
			l.v[i] = -l.v[i]
		} else {
			rand := r.Float64()
			//rand := 0.5
			if de == 4 && rand <= l.lookup[0] {
				l.v[i] = -l.v[i]
			} else if de == 8 && rand <= l.lookup[1] {
				l.v[i] = -l.v[i]
			}
		}
	}
}

func (l *ising) Energy_density() float64 {
	res := 0.
	for i, val := range l.v {
		res -= float64((val * l.v[l.nn[4*i+2]]))
		res -= float64((val * l.v[l.nn[4*i+1]]))
	}
	return res / (float64(l.size * l.size))
}
func (l *ising) Abs_mag_density() float64 {
	res := 0.
	for _, val := range l.v {
		res += float64(val)
	}
	return math.Abs(res / float64(l.size*l.size))
}
