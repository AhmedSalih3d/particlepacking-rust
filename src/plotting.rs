use point;
// PointLabelPosition, Color
use charts::{AxisPosition, Chart, MarkerType, ScaleLinear, ScatterView};

impl charts::PointDatum<f32, f32> for point::Point {
    fn get_x(&self) -> f32 {
        self.0
    }

    fn get_y(&self) -> f32 {
        self.2
    }

    fn get_key(&self) -> String {
        String::new()
    }
}

pub fn plot_points(pini: &[point::Point], name: &str) {
    // Define chart related sizes.
    let width = 800;
    let height = 800;
    let (top, right, bottom, left) = (120, 40, 80, 120);

    // Create a band scale that will interpolate values in [0, 200] to values in the
    // [0, availableWidth] range (the width of the chart without the margins).
    let x = ScaleLinear::new()
        .set_domain(vec![-1.25_f32, 1.25_f32])
        .set_range(vec![0, width - left - right]);

    // Create a linear scale that will interpolate values in [0, 100] range to corresponding
    // values in [availableHeight, 0] range (the height of the chart without the margins).
    // The [availableHeight, 0] range is inverted because SVGs coordinate system's origin is
    // in top left corner, while chart's origin is in bottom left corner, hence we need to invert
    // the range on Y axis for the chart to display as though its origin is at bottom left.
    let y = ScaleLinear::new()
        .set_domain(vec![-1.25_f32, 1.25_f32])
        .set_range(vec![height - top - bottom, 0]);

    // You can use your own iterable as data as long as its items implement the `PointDatum` trait.
    // let scatter_data_1 = vec![(20, 90), (12, 54), (25, 70), (33, 40)];
    // let scatter_data_1 = pini[0 .. 1].to_vec();
    let scatter_data_1 = pini.to_vec();

    // Create Scatter view that is going to represent the data as points.
    let scatter_view_1 = ScatterView::new()
        .set_x_scale(&x)
        .set_y_scale(&y)
        .set_marker_type(MarkerType::Circle)
        // .set_label_position(PointLabelPosition::N)
        // .set_custom_data_label("Apples".to_owned())
        .set_label_visibility(false)
        .load_data(&scatter_data_1)
        .unwrap();

    // Generate and save the chart.
    Chart::new()
        .set_width(width)
        .set_height(height)
        .set_margins(top, right, bottom, left)
        .add_title(String::from("Scatter Chart"))
        .add_view(&scatter_view_1)
        .add_axis_bottom(&x)
        .add_axis_left(&y)
        .add_left_axis_label("Custom Y Axis Label")
        .add_bottom_axis_label("Custom X Axis Label")
        .add_legend_at(AxisPosition::Bottom)
        .save(name)
        .unwrap();
}
