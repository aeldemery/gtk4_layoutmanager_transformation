namespace MathHelper {
    [CCode (cheader_filename = "math_helper.h", cname = "math_helper_perspective_3d")]
    public void perspective_3d (Graphene.Point3D p1,
                                Graphene.Point3D p2,
                                Graphene.Point3D p3,
                                Graphene.Point3D p4,
                                Graphene.Point3D q1,
                                Graphene.Point3D q2,
                                Graphene.Point3D q3,
                                Graphene.Point3D q4,
                                ref Graphene.Matrix m);

    [CCode (cheader_filename = "math_helper.h", cname = "math_helper_singular_value_decomposition")]
    public int singular_value_decomposition (double[] A,
                                             int nrows,
                                             int ncols,
                                             double[] U,
                                             double[] S,
                                             double[] V);

    [CCode (cheader_filename = "math_helper.h", cname = "math_helper_singular_value_decomposition_solve")]
    public void singular_value_decomposition_solve (double[] U,
                                                    double[] S,
                                                    double[] V,
                                                    int nrows,
                                                    int ncols,
                                                    double[] B,
                                                    double[] x);


    public static inline double angle_to_radiants (double angle) {
        return angle * GLib.Math.PI / 180;
    }

    public static inline double map_offset (double x) {
        x = GLib.Math.fmod (x, 180.0);
        if (x < 0.0)
            x += 180.0;

        return x;
    }

    /* Spherical coordinates */
    public static inline double sphere_x (double r, double t, double p) {
        return r * GLib.Math.sin (t) * GLib.Math.cos (p);
    }

    public static inline double sphere_z (double r, double t, double p) {
        return r * GLib.Math.sin (t) * GLib.Math.sin (p);
    }

    public static inline double sphere_y (double r, double t, double p) {
        return r * GLib.Math.cos (t);
    }
}