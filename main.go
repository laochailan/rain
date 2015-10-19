package main

import (
	"fmt"
	"github.com/iand/perlin"
	"math"
	"math/rand"
)

const (
	g  = 9.81
	ρ  = 1000
	γ  = 0.2 / 0.002 // coefficient of friction
	dx = 1e-3

	af          = 0.   // attractive force
	dissipation = 0.05 // coatingheight that is removed from the pane in one second

	ρa = 9e-5 // mean surface energy density of the water coat

	μr = 4e-3 // mean droplet radius
	σr = 2e-3 // std of droplet radius

	μvx   = 2    // mean droplet x velocity
	σv    = 4    // std of droplet velocity
	vterm = 9.65 // terminal droplet velocity

	maxv = 0.4 * μr / 0.002
)

type Droplet struct {
	x, y   float64
	vx, vy float64

	r float64
	h float64

	lifetime float64
}

type Pane struct {
	w, h     int
	fw, fh   float64
	coat     []float64
	dirt     []float64
	droplets []Droplet
}

func lognorm(μ, σ float64) float64 {
	m := math.Log(μ / math.Sqrt(1+σ*σ/μ/μ))
	s := math.Sqrt(math.Log(1 + σ*σ/μ/μ))
	return math.Exp(rand.NormFloat64()*s + m)
}

func (d *Droplet) init() {
	d.x = math.Inf(-1)
	d.y = math.Inf(-1)
	d.r = lognorm(μr, σr)
	d.h = lognorm(d.r, σr)

	d.vx = rand.NormFloat64()*σv + μvx
	μvy := vterm * (1 - math.Exp(-1200*d.r))
	d.vy = rand.NormFloat64()*σv + μvy

	d.lifetime = 0
}

func (d *Droplet) mass() float64 {
	return 2 * d.r * d.r * d.h * ρ
}

func (p *Pane) drawCoating(x, y, r, h float64) (spent float64) {
	ri := int(r + 1)
	xi := int(x + 0.5)
	yi := int(y + 0.5)

	for i := yi - ri; i <= yi+ri; i++ {
		for j := xi - ri; j <= xi+ri; j++ {
			if i < 0 || i >= p.h || j < 0 || j >= p.w {
				continue
			}

			di := float64(i) - y
			dj := float64(j) - x
			if rr := di*di + dj*dj; rr <= r*r {
				vprev := p.coat[i*p.w+j]
				v := 0.1 * h * math.Sqrt(1-rr/r/r)
				if vprev >= v {
					continue
				}
				p.coat[i*p.w+j] = v
				if vprev == 0 {
					spent += dx * dx * dx * h * 0.1
				}
			}
		}
	}

	return spent
}

func (p *Pane) coatForce(x, y, r, h float64) (fx, fy float64) {
	rmax := r + 3
	rmin := 1.0 * r

	ri := int(rmax + 1)
	xi := int(x + 0.5)
	yi := int(y + 0.5)

	c := 0
	for i := yi - ri; i <= yi+ri; i++ {
		for j := xi - ri; j <= xi+ri; j++ {
			if i < 0 || i >= p.h || j < 0 || j >= p.w {
				continue
			}

			di := float64(i) - y
			dj := float64(j) - x
			if rr2 := di*di + dj*dj; rr2 <= rmax*rmax && rr2 > rmin*rmin {
				rr := math.Sqrt(rr2)
				f := 0.0
				if p.coat[i*p.w+j] == 0 {
					f = ρa * p.dirt[i*p.w+j]
				}

				fx += -f * dj / rr / dx
				fy += -f * di / rr / dx
				c++
			}
		}
	}
	if c != 0 {
		fx /= float64(c)
		fy /= float64(c)
	}

	return fx, fy
}

func (p *Pane) update(dt float64) {
	for i := range p.droplets {
		d := &p.droplets[i]
		if d.x+d.r < 0 || d.x-d.r > p.fw || d.y+d.r < -0.4*p.fh || d.y-d.r > p.fh {
			d.init()
			d.x = rand.Float64() * p.fw
			d.y = (1.5*rand.Float64() - 0.4) * p.fh
		}
	}
	for i := range p.droplets {
		d := &p.droplets[i]
		d.lifetime += dt
		fx := 0.0
		fy := d.mass() * g

		cfx, cfy := p.coatForce(d.x/dx, d.y/dx, d.r/dx, d.h/dx)
		fx += cfx
		fy += cfy

		d.vx += dt * fx / d.mass()
		d.vy += dt * fy / d.mass()
		d.vx *= 1 - γ*dt
		d.vy *= 1 - γ*dt

		if v2 := d.vx*d.vx + d.vy*d.vy; v2 > maxv*maxv {
			v := maxv / math.Sqrt(v2)
			d.vx *= v
			d.vy *= v
		}

	}
	for i := range p.droplets {
		d := &p.droplets[i]
		d.x += dt * d.vx
		d.y += dt * d.vy
	}
	for i := range p.droplets {
		d := &p.droplets[i]
		for j := 0; j < i; j++ {
			d2 := &p.droplets[j]
			maxr := d.r
			if d2.r > maxr {
				maxr = d2.r
			}
			ddx := d2.x - d.x
			ddy := d2.y - d.y
			if math.Abs(ddx) > 0.5*maxr || math.Abs(ddy) > 0.5*maxr {
				continue
			}
			if ddx*ddx+ddy*ddy > 0.5*0.5*maxr*maxr {
				continue
			}

			dm := d.mass()
			d2m := d2.mass()
			d.x = (d.x*dm + d2.x*d2m) / (dm + d2m)
			d.y = (d.y*dm + d2.y*d2m) / (dm + d2m)
			maxh := d.h
			if d2.h > maxh {
				maxh = d2.h
			}
			d.r = math.Sqrt(d.r*d.r*d.h/maxh + d2.r*d2.r*d2.h/maxh)
			d.h = maxh

			if d2.vx*d2.vx+d2.vy*d2.vy >= d.vx*d.vx+d.vy*d.vy {
				d.lifetime = d2.lifetime
			}

			d2.init()
		}
	}
	for i, _ := range p.droplets {
		d := &p.droplets[i]
		spent := p.drawCoating(d.x/dx, d.y/dx, d.r/dx*0.7*math.Tanh(d.lifetime), d.h/dx)
		fac := math.Sqrt(d.vx*d.vx + d.vy*d.vy)
		d.h -= 0.04 * fac * spent / d.r / d.r / 2
		if d.h < d.r/10 {
			//oldh := d.h
			d.h = d.r / 10
		}
	}

	for i := range p.coat {
		p.coat[i] -= dissipation * dt * p.dirt[i]
		if p.coat[i] < 0 {
			p.coat[i] = 0
		}
	}

}

func NewPane(w, h, n int) *Pane {
	p := new(Pane)
	p.w = w
	p.h = h
	p.fw = dx * float64(w)
	p.fh = dx * float64(h)
	p.coat = make([]float64, w*h)
	p.dirt = make([]float64, w*h)
	p.droplets = make([]Droplet, n)

	for i := range p.droplets {
		p.droplets[i].init()
	}

	for i := range p.dirt {
		p.dirt[i] = perlin.Noise2D(float64(i%w)/10, float64(i/w)/10, 1, 2, 2, 3)*0.7 + 1
	}
	return p
}

func main() {
	//defer profile.Start(profile.CPUProfile).Stop()
	p := NewPane(1000, 1000, 380)
	fmt.Println("Initialized.")

	for i := 0; i < 1000; i++ {
		for j := 0; j < 200; j++ {
			p.update(0.002)
		}

		err := p.render(fmt.Sprintf("frames/out%05d.png", i))
		fmt.Printf("Frame %d\r", i)
		if err != nil {
			fmt.Println("rendering output:", err)
			return
		}
	}
	fmt.Println("\ndone.")
}
