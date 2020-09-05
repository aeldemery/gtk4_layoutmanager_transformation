public class Gtk4LayoutmanagerTransformation.GlobeTransformationLayout : Gtk.LayoutManager {
    public double position { get; set; }
    public double offset { get; set; }

    /**
     * Those are the 3 Interfaces that we need to implement for our Layout
     *  measure
     *  allocate
     *  get_request_mode
     */

    protected override void measure (Gtk.Widget widget,
                                     Gtk.Orientation orientation,
                                     int for_size,
                                     out int minimum,
                                     out int natural,
                                     out int min_baseline,
                                     out int nat_baseline) {

        int minimum_size = 0;
        int natural_size = 0;

        for (var child = widget.get_first_child ();
             child != null;
             child = child.get_next_sibling ()) {

            int child_min = 0, child_nat = 0;

            if (!child.should_layout ()) {
                continue;
            }

            child.measure (orientation, -1, out child_min, out child_nat, null, null);
            minimum_size = int.max (minimum_size, child_min);
            natural_size = int.max (natural_size, child_nat);
        }

        minimum = minimum_size;
        natural = natural_size * 3;
        min_baseline = -1;
        nat_baseline = -1;
    }

    protected override void allocate (Gtk.Widget widget, int width, int height, int baseline) {
        Gtk.Requisition child_req;
        int j, k;
        float x0, y0;
        float w, h;
        Graphene.Point3D p1, p2, p3, p4;
        Graphene.Point3D q1, q2, q3, q4;
        double t_1, t_2, p_1, p_2;
        double r;
        Graphene.Matrix m = {};
        Gsk.Transform transform = new Gsk.Transform ();

        /* for simplicity, assume all children are the same size */
        widget.get_preferred_size (out child_req, null);
        w = child_req.width;
        h = child_req.height;

        r = 300;
        x0 = y0 = 300;

        for (var child = widget.get_first_child (), i = 0;
             child != null;
             child = child.get_next_sibling (), i++) {
            j = i / 36;
            k = i % 36;

            child.set_child_visible (false);

            p1 = { w, h, 1f };
            p2 = { w, 0f, 1f };
            p3 = { 0f, 0f, 1f };
            p4 = { 0f, h, 1f };

            t_1 = MathHelper.angle_to_radiants (MathHelper.map_offset (offset + 10 * j));
            t_2 = MathHelper.angle_to_radiants (MathHelper.map_offset (offset + 10 * (j + 1)));
            p_1 = MathHelper.angle_to_radiants (position + 10 * k);
            p_2 = MathHelper.angle_to_radiants (position + 10 * (k + 1));

            if (t_2 < t_1) {
                continue;
            }

            if (MathHelper.sphere_z (r, t_1, p_1) > 0 ||
                MathHelper.sphere_z (r, t_2, p_1) > 0 ||
                MathHelper.sphere_z (r, t_1, p_2) > 0 ||
                MathHelper.sphere_z (r, t_2, p_2) > 0) {
                continue;
            }

            child.set_child_visible (true);

            q1 = {
                (float) (x0 + MathHelper.sphere_x (r, t_1, p_1)),
                (float) (y0 + MathHelper.sphere_y (r, t_1, p_1)),
                (float) MathHelper.sphere_z (r, t_1, p_1)
            };

            q2 = {
                (float) (x0 + MathHelper.sphere_x (r, t_2, p_1)),
                (float) (y0 + MathHelper.sphere_y (r, t_2, p_1)),
                (float) MathHelper.sphere_z (r, t_2, p_1)
            };

            q3 = {
                (float) (x0 + MathHelper.sphere_x (r, t_2, p_2)),
                (float) (y0 + MathHelper.sphere_y (r, t_2, p_2)),
                (float) MathHelper.sphere_z (r, t_2, p_2)
            };

            q4 = {
                (float) (x0 + MathHelper.sphere_x (r, t_1, p_2)),
                (float) (y0 + MathHelper.sphere_y (r, t_1, p_2)),
                (float) MathHelper.sphere_z (r, t_1, p_2)
            };

            /* Get a matrix that moves p1 -> q1, p2 -> q2, ... */
            MathHelper.perspective_3d (p1, p2, p3, p4,
                                       q1, q2, q3, q4,
                                       ref m);

            transform = transform.matrix (m);

            /* Since our matrix was built for transforming points with z = 1,
             * prepend a translation to the z = 1 plane.
             */
            transform = transform.translate_3d ({ 0, 0, 1 });

            child.allocate ((int) w, (int) h, -1, transform);
        }
    }

    protected override Gtk.SizeRequestMode get_request_mode (Gtk.Widget widget) {
        return Gtk.SizeRequestMode.CONSTANT_SIZE;
    }
}