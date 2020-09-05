public class Gtk4LayoutmanagerTransformation.GlobeWidget : Gtk.Widget {
    int64 start_time;
    int64 end_time;

    double start_position;
    double end_position;

    double start_offset;
    double end_offset;

    bool animating;

    public GlobeWidget () {
        this.install_action ("rotate", "(ii)", (Gtk.WidgetActionActivateFunc)rotate_sphere);
        this.add_binding_action (Gdk.Key.Left, 0, "rotate", "(ii)", Gtk.Orientation.HORIZONTAL, -1);
        this.add_binding_action (Gdk.Key.Right, 0, "rotate", "(ii)", Gtk.Orientation.HORIZONTAL, 1);
        this.add_binding_action (Gdk.Key.Up, 0, "rotate", "(ii)", Gtk.Orientation.VERTICAL, 1);
        this.add_binding_action (Gdk.Key.Down, 0, "rotate", "(ii)", Gtk.Orientation.VERTICAL, -1);

        /* here is where we use our custom layout manager */
        this.set_layout_manager_type (typeof (GlobeTransformationLayout));
        this.focusable = true;
    }

    public void add_child (Gtk.Widget child) {
        child.set_parent (this);
    }

    /* From clutter-easing.c, based on Robert Penner's
     * infamous easing equations, MIT license.
     */
    double ease_out_cubic (double t) {
        double p = t - 1;

        return p * p * p + 1;
    }

    bool update_position (Gtk.Widget unused, Gdk.FrameClock clock) {
        var layout = (GlobeTransformationLayout) this.layout_manager;

        var now = clock.get_frame_time ();

        if (now >= this.end_time) {
            this.animating = false;
            return GLib.Source.REMOVE;
        }

        double t = (now - this.start_time) / (this.end_time - this.start_time);

        t = ease_out_cubic (t);

        layout.position = this.start_position + t * (this.end_position - this.start_position);
        layout.offset = this.start_offset + t * (this.end_offset - this.start_offset);

        this.queue_allocate ();

        return GLib.Source.CONTINUE;
    }

    void rotate_sphere (Gtk.Widget self, string action, GLib.Variant parameters) {
        var layout = (GlobeTransformationLayout) this.layout_manager;
        Gtk.Orientation orientation;
        int direction;

        parameters.get ("(ii)", out orientation, out direction);

        this.end_position = this.start_position = layout.position;
        this.end_offset = this.start_offset = layout.offset;

        if (orientation == Gtk.Orientation.HORIZONTAL) {
            this.end_position += 10 * direction;
        } else {
            this.end_offset += 10 * direction;
        }
        this.start_time = GLib.get_monotonic_time ();
        this.end_time = (int64) (this.start_time + 0.5 * GLib.TimeSpan.SECOND);

        if (!this.animating) {
            this.add_tick_callback (update_position);
            this.animating = true;
        }
    }

    protected override void snapshot (Gtk.Snapshot snapshot) {
        for (var child = this.get_first_child ();
             child != null;
             child = child.get_next_sibling ()) {
            /* our layout manager sets this for children that are out of view */
            if (!child.get_child_visible ()) {
                continue;
            }

            this.snapshot_child (child, snapshot);
        }
    }

    protected override void dispose () {
        var child = this.get_first_child ();
        while (child != null) {
            child.unparent ();
        }
        base.dispose ();
    }
}
