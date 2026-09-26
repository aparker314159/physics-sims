
use std::error::Error;

use cairo::{Context, ImageSurface};
use plotters::coord::Shift;
use plotters::prelude::*;
use plotters_cairo::CairoBackend;

pub const WIDTH: usize = 1200;
pub const HEIGHT: usize = 800;


pub struct Graph2d<'w> {
    x_bounds: (f64, f64),
    y_bounds: (f64, f64),
    surface: &'w mut ImageSurface,
    root: DrawingArea<CairoBackend<'w>, Shift>,
    pixels: Vec<u32>,
}


impl<'w> Graph2d<'w> {
    pub fn new(
        x_bounds: (f64,f64), y_bounds: (f64,f64), 
        surface: &'w mut ImageSurface, context: &'w Context) 
        -> Result<Self, Box<dyn Error>> {
        let root = CairoBackend::new(&context, (WIDTH as u32, HEIGHT as u32))?.into_drawing_area();

        Ok(Graph2d { 
            x_bounds, 
            y_bounds, 
            surface, 
            root, 
            pixels: vec![0; WIDTH * HEIGHT] 
        })
    }

    pub fn draw_frame(&mut self, data: Vec<Vec<(f64, f64)>>) -> Result<(), Box<dyn Error>> {
        self.root.fill(&WHITE)?;
        let mut chart = ChartBuilder::on(&self.root)
            .margin(20)
            .set_all_label_area_size(50)
            .build_cartesian_2d(self.x_bounds.0 .. self.x_bounds.1, 
                                self.y_bounds.0 .. self.y_bounds.1)?;

        chart.configure_mesh()
            .light_line_style(RGBColor(127, 127, 127).stroke_width(1))
            .bold_line_style(BLACK.stroke_width(1))
            .x_desc("x")
            .y_desc("y")
            .x_labels(20)
            .y_labels(10)
            .axis_style(BLACK.stroke_width(4))
            .draw()?;

        let colors = vec![BLUE.stroke_width(3), RED.stroke_width(3), GREEN.stroke_width(3), YELLOW.stroke_width(3)];
        let mut color_it = colors.iter().cycle();
        for series in data {
            chart.draw_series(LineSeries::new(
                series,
                color_it.next().expect("").stroke_width(3),
            ))?;
        }

        self.root.present()?;

        let stride = self.surface.stride() as usize;
        self.surface.with_data(|data| {
            for (row, output) in data.chunks_exact(stride).zip(self.pixels.chunks_exact_mut(WIDTH)) {
                for (bytes, pixel) in row[..WIDTH * 4].chunks_exact(4).zip(output) {
                    *pixel = u32::from_ne_bytes([bytes[0], bytes[1], bytes[2], bytes[3]]);
                }
            }
        })?;

        return Ok(());
    }

    pub fn get_pixels(&self) -> &[u32] {
        self.pixels.as_slice()
    }
}

