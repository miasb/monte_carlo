package markov

import (
	"math"
	"math/rand"
)

type Markov struct {
	Energy []float64
	Mag    []float64
	mv     int
}

func New(len int, eq float64) *Markov {
	chains := Markov{}

	chains.mv = int(float64(len) * eq / 100.)
	chains.Energy = make([]float64, len)
	chains.Mag = make([]float64, len)

	return &chains
}

func (chains *Markov) Av_Energy() float64 {
	return avg(chains.Energy[chains.mv:])
}
func (chains *Markov) Av_Magnetization() float64 {
	return avg(chains.Mag[chains.mv:])
}
func (chains *Markov) Magn_Suscep(t float64, size int) float64 {
	dev := std_dev(chains.Mag[chains.mv:])

	return 1 / t * dev * dev * float64(size*size)
}
func (chains *Markov) Spec_Heat(t float64, size int) float64 {
	dev := std_dev(chains.Energy[chains.mv:])

	return 1 / (t * t) * dev * dev * float64(size*size)
}

func avg(v []float64) float64 {
	res := 0.

	for _, num := range v {
		res += float64(num)
	}
	return res / float64(len(v))
}
func std_dev(v []float64) float64 {
	res := 0.
	avg := avg(v)

	for _, num := range v {
		res += math.Pow(num-avg, 2)
	}
	return math.Sqrt(res / float64(len(v)-1))
}

func (chains *Markov) Std_Err_Primary(chain []float64) float64 {
	variance := std_dev(chain[chains.mv:])
	length := float64(len(chain) - chains.mv)

	return math.Sqrt(2*chains.Autocorr_Time(chain)/length) * variance
}

func (chains *Markov) Autocorr_Time(chain []float64) float64 {
	eq_len := len(chains.Energy) - chains.mv

	dev := std_dev(chain[chains.mv:])

	C0 := dev * dev
	Ct := 0.

	tau := 0.
	rho := 0.
	t := 1

	for rho >= 0 {
		tau += rho

		Ct = 0
		ymin := avg(chain[chains.mv:])
		yplus := avg(chain[chains.mv+t:])

		for i := 0; i < eq_len-t; i++ {
			Ct += (chain[i+chains.mv] - ymin) * (chain[i+chains.mv+t] - yplus)
		}

		rho = Ct / float64(eq_len-t) / C0
		t += 1
	}
	return tau + 0.5
}

func (chains *Markov) Error_Spec_Heat(algo int, size int, t float64) float64 {
	nb := 20
	M := 1000
	switch algo {
	case 1:
		return chains.blocking_c(nb, size, t)
	case 2:
		return chains.bootstrap_c(M, size, t)
	default:
		return chains.error_prop_c(size, t)
	}
}
func (chains *Markov) error_prop_c(size int, t float64) float64 {
	av := chains.Av_Energy()
	length := len(chains.Energy) - chains.mv
	err_prop_chain := make([]float64, len(chains.Energy))

	for i := 0; i < length; i++ {
		X := chains.Energy[i+chains.mv]
		err_prop_chain[i+chains.mv] = X*X - 2*av*X
	}

	return chains.Std_Err_Primary(err_prop_chain) * float64(size*size) / (t * t)
}
func (chains *Markov) blocking_c(nb int, size int, t float64) float64 {
	length := len(chains.Energy) - chains.mv

	blocksize := length / nb
	block_meas := make([]float64, nb)

	for i := 0; i < nb; i++ {
		dev := std_dev(chains.Energy[i*blocksize+chains.mv : (i+1)*blocksize+chains.mv])
		block_meas[i] = 1 / (t * t) * dev * dev * float64(size*size)
	}

	return std_dev(block_meas[1:]) / math.Sqrt(float64(nb)-1.)
}
func (chains *Markov) bootstrap_c(M int, size int, t float64) float64 {
	length := len(chains.Energy) - chains.mv

	tau := chains.Autocorr_Time(chains.Energy)
	Nindep := int(float64(len(chains.Energy)) / (2 * tau))

	pseudo_chain := make([]float64, M)

	for i := 0; i < M; i++ {
		bootstrap := make([]float64, Nindep)
		for j := 0; j < Nindep; j++ {
			rand := rand.Intn(length)
			bootstrap[j] = chains.Energy[rand+chains.mv]
		}
		dev := std_dev(bootstrap)
		pseudo_chain[i] = 1. / (t * t) * dev * dev * float64(size*size)
	}
	return std_dev(pseudo_chain)
}

func (chains *Markov) Error_Magn_Suscep(algo int, size int, t float64) float64 {
	nb := 20
	M := 1000
	switch algo {
	case 1:
		return chains.blocking_x(nb, size, t)
	case 2:
		return chains.bootstrap_x(M, size, t)
	default:
		return chains.error_prop_x(size, t)
	}
}
func (chains *Markov) error_prop_x(size int, t float64) float64 {
	av := chains.Av_Magnetization()
	length := len(chains.Mag) - chains.mv
	err_prop_chain := make([]float64, len(chains.Mag))

	for i := 0; i < length; i++ {
		X := chains.Mag[i+chains.mv]
		err_prop_chain[i+chains.mv] = X*X - 2*av*X
	}

	return chains.Std_Err_Primary(err_prop_chain) * float64(size*size) / t
}
func (chains *Markov) blocking_x(nb int, size int, t float64) float64 {
	length := len(chains.Mag) - chains.mv

	blocksize := length / nb
	block_meas := make([]float64, nb)

	for i := 0; i < nb; i++ {
		dev := std_dev(chains.Mag[i*blocksize+chains.mv : (i+1)*blocksize+chains.mv])
		block_meas[i] = 1 / t * dev * dev * float64(size*size)
	}

	return std_dev(block_meas[1:]) / math.Sqrt(float64(nb)-1)
}
func (chains *Markov) bootstrap_x(M int, size int, t float64) float64 {
	length := len(chains.Mag) - chains.mv

	tau := chains.Autocorr_Time(chains.Mag)
	Nindep := int(float64(len(chains.Mag)) / (2 * tau))

	pseudo_chain := make([]float64, M)

	for i := 0; i < M; i++ {
		bootstrap := make([]float64, Nindep)
		for j := 0; j < Nindep; j++ {
			rand := rand.Intn(length)
			bootstrap[j] = chains.Mag[rand+chains.mv]
		}
		dev := std_dev(bootstrap)
		pseudo_chain[i] = 1. / t * dev * dev * float64(size*size)
	}
	return std_dev(pseudo_chain)
}
