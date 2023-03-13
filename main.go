package main

import (
	"fmt"
	"math/rand"
	"monte_carlo/markov"
	"monte_carlo/xy"
	"os"
)

func check(e error) {
	if e != nil {
		panic(e)
	}
}

func main() {
	temp := 0.9
	size := 16
	length := 60000
	n := 2
	dt := 0.046

	l := xy.New(size, temp, n, dt)
	chains := markov.New(length, 20)
	r := rand.New(rand.NewSource(99))
	write_to_file := false

	for temp < 3.0 {

		folder := "../data/indep/high_t/"
		path := folder + fmt.Sprintf("%.3f", dt) + ".dat"
		f, err := os.Create(path)
		check(err)
		defer f.Close()

		for i := range chains.Energy {
			chains.Energy[i] = l.Energy_density()
			chains.Mag[i] = l.Abs_mag_density()
			l.Sweep_HMC(r)
		}
		//av_energy := chains.Av_Energy()
		av_mag := chains.Av_Magnetization()
		//mag_sus := chains.Magn_Suscep(temp, size)
		//spec_h := chains.Spec_Heat(temp, size)
		if write_to_file {
			string := fmt.Sprintf("%d", n) + " " +
				fmt.Sprintf("%f", 2*chains.Autocorr_Time(chains.Energy)*float64(n)) + "\n"
			f.WriteString(string)
		} else {
			fmt.Println(l.Get_t(), av_mag, chains.Std_Err_Primary(chains.Mag))
		}
		temp += 0.05
		l.Set_t(temp)
		l.Cold_Start()

		if write_to_file {
			f.Sync()
		}
	}

	// l := ising.New(size, temp)
	// chains := markov.New(length, 10)
	// r := rand.New(rand.NewSource(99))

	// for j := 0; j < 1; j++ {
	// 	for i := 0; i < length; i++ {
	// 		chains.Energy[i] = l.Energy_density()
	// 		chains.Mag[i] = l.Abs_mag_density()
	// 		l.Sweep(r)
	// 	}
	// 	av_energy := chains.Av_Energy()
	// 	av_mag := chains.Av_Magnetization()
	// 	mag_sus := chains.Magn_Suscep(temp, size)
	// 	spec_h := chains.Spec_Heat(temp, size)
	// 	fmt.Println(temp, av_energy, av_mag, mag_sus, spec_h,
	// 		chains.Std_Err_Primary(chains.Energy),
	// 		chains.Std_Err_Primary(chains.Mag),
	// 		chains.Error_Spec_Heat(0, size, temp),
	// 		chains.Error_Magn_Suscep(0, size, temp))
	//
	// 	temp += 0.1
	// 	l.Set_t(temp)
	// }
}
