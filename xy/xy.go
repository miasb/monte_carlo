package xy

import (
	"math"
	"math/rand"
)

type xy struct {
	size   int
	t      float64
	n      int
	dt     float64
	v      []float64
	v_bck  []float64
	nn     []int
	lookup []float64
}

func New(size int, t float64, n int, dt float64) *xy {
	l := xy{size: size, t: t, n: n, dt: dt}
	l.v = make([]float64, size*size)
	l.v_bck = make([]float64, size*size)
	for i := range l.v {
		l.v[i] = 1.
	}

	l.nn = make([]int, 4*size*size)
	calc_nn(l.nn, size)

	l.lookup = make([]float64, 16)
	for i := range l.lookup {
		l.lookup[i] = float64(i) * 2. * math.Pi / 16.
	}

	return &l
}

func (l *xy) Set_t(temp float64) {
	l.t = temp
}
func (l *xy) Get_t() float64 {
	return l.t
}
func (l *xy) Set_n(n int) {
	l.n = n
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

func (l *xy) dE(random float64, i int) float64 {
	res := math.Cos(random - l.v[l.nn[4*i+2]])
	res += math.Cos(random - l.v[l.nn[4*i+1]])
	res += math.Cos(random - l.v[l.nn[4*i+3]])
	res += math.Cos(random - l.v[l.nn[4*i]])
	return res
}

func (l *xy) Sweep(r *rand.Rand) {
	for i := range l.v {
		random := r.Float64() * 2 * math.Pi
		de := l.dE(l.v[i], i) - l.dE(random, i)
		if de <= 0 {
			l.v[i] = random
		} else {
			r_acceptance := r.Float64()
			if r_acceptance <= math.Exp(-de/l.t) {
				l.v[i] = random
			}
		}
	}
}

func (l *xy) Sweep_HMC(r *rand.Rand) {
	// the new config is accepted with a probability
	// so the old one is stored here
	copy(l.v_bck, l.v)

	H_start := l.calc_H()
	l.calc_trajectory(r)
	H_end := l.calc_H()
	//fmt.Println(math.Exp((-H_end + H_start) / l.t))

	random := r.Float64()
	if random > math.Exp(-(H_end-H_start)/l.t) {
		copy(l.v, l.v_bck)
	}
}

func (l *xy) calc_H() float64 {
	sum := 0.

	for i, val := range l.v {
		sum += math.Cos(val - l.v[l.nn[4*i+2]])
		sum += math.Cos(val - l.v[l.nn[4*i+1]])
	}

	return -sum
}

func log_normal(x float64) float64 {
	//return 1 / (x * math.Sqrt(2*math.Pi)) *
	//	math.Exp(
	//		-(math.Pow(math.Log(x), 2))/2)
	return 1 / math.Sqrt(2*math.Pi) *
		math.Exp(-0.5*x*x)
}
func (l *xy) grad(x float64) float64 {
	fx := log_normal(x)
	fxh := log_normal(x + l.dt)

	return (fxh - fx) / l.dt
}

func (l *xy) calc_trajectory(r *rand.Rand) {
	momenta := make([]float64, l.size*l.size)

	for i := 0; i < l.size*l.size; i++ {
		momenta[i] = r.NormFloat64()
		// initial half-step
		// Calculate the gradient of the normal distribution
		Vq := l.grad(l.v[i])
		momenta[i] = momenta[i] - Vq*l.dt/2
	}

	for i := 0; i < l.n; i++ {
		for j := 0; j < l.size*l.size; j++ {
			l.v[j] += momenta[j] * l.dt
			Vq := l.grad(l.v[j])
			momenta[j] -= Vq * l.dt
		}
	}
}

func (l *xy) Abs_mag_density() float64 {
	cos_part := 0.0
	sin_part := 0.0

	for _, val := range l.v {
		cos_part += math.Cos(val)
		sin_part += math.Sin(val)
	}

	return math.Sqrt(cos_part*cos_part+sin_part*sin_part) / (float64(l.size * l.size))
}
func (l *xy) Energy_density() float64 {
	//fmt.Println(l.v[123])
	return l.calc_H() / float64(l.size*l.size)
}

func (l *xy) Cold_Start() {
	for i := range l.v {
		l.v[i] = 0.
	}
}
