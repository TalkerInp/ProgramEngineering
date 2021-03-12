package main

//Variant 6

import (
	"fmt"
	"math"
	"strconv"

	"gonum.org/v1/gonum/mat"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/plotutil"
	"gonum.org/v1/plot/vg"
)

type defusion struct {
	D   float64
	t   float64
	Re  float64
	Pr  float64
	cs  float64
	cX0 float64
}

type item struct {
	name  string
	value plot.Thumbnailer
}

type pipeline struct {
	l          float64
	d          float64
	mydefusion defusion
}

func matPrint(X mat.Matrix) {
	fa := mat.Formatted(X, mat.Prefix(""), mat.Squeeze())
	fmt.Printf("%v\n", fa)
}

func liab2Solution(mypipe pipeline) *mat.Dense {
	dx := 0.01
	var dt float64
	dt = 50
	XN := int((mypipe.l-0)/dx + 1)
	TN := int((mypipe.mydefusion.t-0)/dt + 1)

	C := mat.NewDense(XN, TN, nil) // объявление-обнуление матрицы

	for i := 0; i < XN; i++ {
		answer := 200 + 50*float64(i)*dx
		C.Set(i, 0, answer)
	}

	for i := 0; i < TN; i++ {
		C.Set(0, i, 200)
	}

	a := -(mypipe.mydefusion.D * dt) / math.Pow(dx, 2)
	b := 1 + 2*((mypipe.mydefusion.D*dt)/math.Pow(dx, 2))
	c := a
	nu := mypipe.mydefusion.D * (1 + 0.5*(0.55*math.Pow(mypipe.mydefusion.Re, 0.5)*math.Pow(mypipe.mydefusion.Pr, (1/3)))) / mypipe.d

	for j := 0; j < TN-1; j++ {
		alfa := make([]float64, 0, XN)
		alfa = append(alfa, 0)
		beta := make([]float64, 0, XN)
		beta = append(beta, 200)
		for z := 1; z <= XN-1; z++ {
			alfa = append(alfa, -(a / (b + c*alfa[z-1])))
			beta = append(beta, (C.At(z, j)-c*beta[z-1])/(b+c*alfa[z-1]))
		}
		solution := (nu*mypipe.mydefusion.cs*dx + beta[XN-2]) / (1 - alfa[XN-2] + nu*dx)
		C.Set(XN-1, j+1, solution)
		for m := XN - 2; m >= 0; m-- {
			answer := alfa[m]*C.At(m+1, j+1) + beta[m]
			C.Set(m, j+1, answer)
		}
	}
	matPrint(C)
	return C
}

func main() {
	var mypipe = pipeline{
		l: 0.2,
		d: 0.01,
		mydefusion: defusion{
			D:   3e-6,
			t:   10000,
			Re:  6,
			Pr:  400,
			cs:  10,
			cX0: 200,
		},
	}
	C := liab2Solution(mypipe)
	makePlot(C, mypipe)

}

func makePlot(C *mat.Dense, mypipe pipeline) {
	p, err := plot.New()
	if err != nil {
		panic(err)
	}
	dx := 0.01
	rows, columns := C.Caps()
	p.Title.Text = "Laba2"
	p.X.Label.Text = "X"
	p.Y.Label.Text = "Y"

	var ps []plot.Plotter
	var items []item
	name := ""
	pts := make(plotter.XYs, rows)
	for j := 0; j < columns; j++ {
		for i := 0; i < rows; i++ {
			pts[i].X = mypipe.l + dx*float64(i)
			pts[i].Y = C.At(i, j)
		}
		name = strconv.Itoa(j)
		l, err := plotter.NewLine(pts)
		if err != nil {
			fmt.Println("Error")
		}
		l.Color = plotutil.Color(j)
		ps = append(ps, l)
		if name != "" {
			items = append(items, item{name: name, value: l})
			name = ""
		}
		p.Add(ps...)
		for _, v := range items {
			p.Legend.Add(v.name, v.value)
		}
	}

	if err := p.Save(15*vg.Inch, 15*vg.Inch, "pointSLaba2.png"); err != nil {
		panic(err)
	}

}
