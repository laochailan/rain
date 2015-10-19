package main

import (
	"image"
	"image/color"
	"image/png"
	"math"
	"os"
)

func drawDroplet(img *image.Gray, x, y, r, h float64) {
	ri := int(r + 1)
	xi := int(x + 0.5)
	yi := int(y + 0.5)
	for i := yi - ri; i <= yi+ri; i++ {
		for j := xi - ri; j <= xi+ri; j++ {
			if i < 0 || i >= img.Rect.Dy() || j < 0 || j >= img.Rect.Dx() {
				continue
			}

			di := float64(i) - y
			dj := float64(j) - x
			if rr := di*di + dj*dj; rr <= r*r {
				vprev := img.GrayAt(j, i).Y
				v := 100*h/r*math.Pow(1-rr/r/r, 1.0)/math.Sqrt(1-di/r) + float64(vprev)
				if v > 255 {
					v = 255
				}

				img.SetGray(j, i, color.Gray{uint8(v)})
			}
		}
	}
}

func (p *Pane) render(filename string) error {
	w := p.w
	h := p.h
	mtp := 1 / dx
	img := image.NewGray(image.Rect(0, 0, w, h))

	for y := 0; y < h; y++ {
		for x := 0; x < w; x++ {

			img.SetGray(x, y, color.Gray{uint8(math.Sqrt(p.coat[y*w+x]) * 80)})
		}
	}

	for _, d := range p.droplets {
		drawDroplet(img, d.x*mtp, d.y*mtp, d.r*mtp, d.h*mtp)
	}

	file, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	err = png.Encode(file, img)
	if err != nil {
		return err
	}

	return nil
}
