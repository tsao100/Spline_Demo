#include <X11/Xlib.h>
#include "tinyspline.h"
#include <stdlib.h>
#include <stdio.h>
#include <X11/keysym.h>

Display *display;
Window window;
GC gc;

// Spline instance
static tsBSpline spline;
static int spline_ready = 0;

// 初始化視窗
void ps_init_window_() {
    display = XOpenDisplay(NULL);
    Window root = DefaultRootWindow(display);
    window = XCreateSimpleWindow(display, root, 50,50,800,600,1,0x000000,0xFFFFFF);
    XSelectInput(display, window, ExposureMask | ButtonPressMask);
    XMapWindow(display, window);
    gc = XCreateGC(display, window, 0, NULL);
    XFlush(display);
}

// 畫點
void ps_draw_point_(int *x, int *y) {
    XColor red_color, exact;
    Colormap cmap = DefaultColormap(display, 0);

    if (XAllocNamedColor(display, cmap, "red", &red_color, &exact)) {
        XSetForeground(display, gc, red_color.pixel);

        XDrawArc(display, window, gc, *x - 3, *y - 3, 6, 6, 0, 360 * 64);
        XFlush(display);

        // 你也可以在這裡恢復原來的顏色，例如白色或黑色：
        XSetForeground(display, gc, BlackPixel(display, 0));
    } else {
        fprintf(stderr, "Failed to allocate color 'red'\n");
    }
}




// 畫 Catmull-Rom spline
void ps_draw_spline_(int *x, int *y, int *n)
{
    int i;
    int dim = 2;
    tsStatus status;
    tsBSpline curve;
    tsError err;

    // 將 Fortran 傳來的 (x, y) 組成 points[]
    size_t n_points = (size_t)*n;
    tsReal *points = malloc(sizeof(tsReal) * n_points * dim);
    if (!points) return;

    XClearWindow(display, window);

    for (i = 0; i < *n; i++) {
        points[i * 2 + 0] = (tsReal)x[i];
        points[i * 2 + 1] = (tsReal)y[i];
        ps_draw_point_(&x[i], &y[i]);
//        printf("(%d %d).\n",x[i], y[i]);
    }

    // 產生 Catmull-Rom spline
    err = ts_bspline_interpolate_catmull_rom(points, n_points, dim,
                                             0.5, NULL, NULL, 1e-4,
                                             &curve, &status);
    free(points);

    if (err != TS_SUCCESS) {
        fprintf(stderr, "Spline interpolation failed: %s\n", status.message);
        return;
    }

    // 取樣 spline
    tsReal *samples = NULL;
    size_t sample_count = 100;
    size_t actual_count;

    err = ts_bspline_sample(&curve, sample_count, &samples, &actual_count, &status);
    if (err != TS_SUCCESS) {
        fprintf(stderr, "Spline sampling failed: %s\n", status.message);
        ts_bspline_free(&curve);
        return;
    }

    // 繪圖
    for (i = 0; i < sample_count - 1; i++) {
        int x1 = (int)samples[i * 2 + 0];
        int y1 = (int)samples[i * 2 + 1];
        int x2 = (int)samples[(i + 1) * 2 + 0];
        int y2 = (int)samples[(i + 1) * 2 + 1];
        XDrawLine(display, window, gc, x1, y1, x2, y2);
    }
    XFlush(display);
    // 清理
    free(samples);
    ts_bspline_free(&curve);
}

// 等待滑鼠點擊
void ps_wait_click_(int *x, int *y) {
    XEvent ev;
    while (1) {
        XNextEvent(display, &ev);
        if (ev.type == ButtonPress) {
            if (ev.xbutton.button == Button3) {  // Right-click to quit
                *x = -1;
                *y = -1;
                break;
            } else {
            	if (ev.xbutton.button == Button1) {
    		   *x = ev.xbutton.x;
                   *y = ev.xbutton.y;
		}
                break;
            }
        } else if (ev.type == KeyPress) {
            KeySym keysym = XLookupKeysym(&ev.xkey, 0);
            if (keysym == XK_Return) {  // Enter key
                *x = -1;
                *y = -1;
                break;
            }
        }
    }
}


// 建立 spline
void ps_spline_new_(int *deg, int *dim, int *n_ctrlp, int *type) {
    tsStatus status;
    tsError err = ts_bspline_new((size_t)*n_ctrlp,
                                 (size_t)*dim,
                                 (size_t)*deg,
                                 (tsBSplineType)*type,
                                 &spline,
                                 &status);
    spline_ready = (err == TS_SUCCESS);
    if (!spline_ready) {
        fprintf(stderr, "Error: %s\n", status.message);
    }
}

// 複製 spline
void ps_spline_copy_(tsBSpline *dest) {
    if (!spline_ready) return;
    ts_bspline_copy(&spline, dest, NULL);
}

// 釋放 spline
void ps_spline_free_() {
    if (spline_ready) ts_bspline_free(&spline);
    spline_ready = 0;
}

// 插入控制點
void ps_spline_set_ctrlp_(double *ctrlp) {
    if (!spline_ready) return;
    ts_bspline_set_control_points(&spline, ctrlp, NULL);
}

// 讀取控制點
void ps_spline_get_ctrlp_(double *ctrlp_out) {
    if (!spline_ready) return;
    tsReal *cp = NULL;
    ts_bspline_control_points(&spline, &cp, NULL);
    size_t n = ts_bspline_len_control_points(&spline);
    for (size_t i = 0; i < n; i++) ctrlp_out[i] = cp[i];
}

// 使用預設 knots 自動設
void ps_spline_set_knots_(double *knots_in, int *n_knots, int *status_out) {
    if (!spline_ready) return;

    tsStatus status;
    tsError err = ts_bspline_set_knots(&spline,
                                       (const tsReal *)knots_in,
                                       &status);
    if (status_out) *status_out = (err == TS_SUCCESS) ? 0 : 1;

    if (err != TS_SUCCESS)
        fprintf(stderr, "Set knots failed: %s\n", status.message);
}

// 插入 knot
void ps_spline_insert_knot_(double *u, int *n, int *k_out) {
    if (!spline_ready) return;
    size_t k;
    ts_bspline_insert_knot(&spline, *u, (size_t)*n, &spline, &k, NULL);
    *k_out = (int)k;
}

// 拆解為 Bezier
void ps_spline_to_beziers_(tsBSpline *out) {
    if (!spline_ready) return;
    ts_bspline_to_beziers(&spline, out, NULL);
}

// Spline 插值
void ps_spline_interpolate_cubic_natural_(double *points, int *n, int *dim, int *status_out) {
    tsStatus status;
    tsError err = ts_bspline_interpolate_cubic_natural(
        (const tsReal *)points,
        (size_t)*n,
        (size_t)*dim,
        &spline,
        &status
    );

    spline_ready = (err == TS_SUCCESS);
    if (status_out)
        *status_out = (err == TS_SUCCESS) ? 0 : 1;

    if (err != TS_SUCCESS)
        fprintf(stderr, "Interpolation failed: %s\n", status.message);
}

void ps_spline_interpolate_catmull_rom_(double *points, int *n, int *dim,
                                    double *alpha, double *epsilon,
                                    int *status_out)
{
    tsStatus status;
    tsError err;

    // 你可以不提供 first/last，即可建立「非接觸」端點的自然 Catmull-Rom spline
    err = ts_bspline_interpolate_catmull_rom(
        (const tsReal *)points,
        (size_t)*n,
        (size_t)*dim,
        (tsReal)*alpha,
        NULL,   // first
        NULL,   // last
        (tsReal)*epsilon,
        &spline,
        &status
    );

    spline_ready = (err == TS_SUCCESS);
    if (status_out)
        *status_out = (err == TS_SUCCESS) ? 0 : 1;

    if (err != TS_SUCCESS)
        fprintf(stderr, "Catmull-Rom interpolation failed: %s\n", status.message);
}


// 評估單一 u
void ps_spline_eval_(double *u, double *out) {
    if (!spline_ready) return;
    tsDeBoorNet net = ts_deboornet_init();
    ts_bspline_eval(&spline, *u, &net, NULL);
    size_t d = ts_deboornet_dimension(&net);
    const tsReal *pts = ts_deboornet_points_ptr(&net);
    for (size_t i=0;i<d;i++) out[i] = pts[i];
    ts_deboornet_free(&net);
}

// 取樣多點
void ps_spline_sample_(int *samples, double *out, int *actual) {
    if (!spline_ready) return;
    size_t act = 0;
    ts_bspline_sample(&spline, (size_t)*samples, &out, &act, NULL);
    *actual = (int)act;
}

// 分割 splines
void ps_spline_split_(double *u, tsBSpline *out, int *k_out) {
    if (!spline_ready) return;
    size_t k;
    ts_bspline_split(&spline, *u, out, &k, NULL);
    *k_out = (int)k;
}

// 計算字元數量
void ps_spline_degree_(int *deg_out) {
    if (!spline_ready) return;
    *deg_out = (int)ts_bspline_degree(&spline);
}

// 計算 knot range
void ps_spline_domain_(double *min, double *max) {
    if (!spline_ready) return;
    ts_bspline_domain(&spline, min, max);
}

// 是否 Closed
void ps_spline_is_closed_(double *eps, int *closed_out) {
    if (!spline_ready) return;
    ts_bspline_is_closed(&spline, *eps, closed_out, NULL);
}

// 計算 chord lengths
void ps_spline_chord_lengths_(double *knots, int *num, double *lengths) {
    ts_bspline_chord_lengths(&spline, knots, (size_t)*num, lengths, NULL);
}

// 等距 knot seq
void ps_spline_equidistant_knots_(double *knots, int *num_knot_seq, int *num_samples) {
    ts_chord_lengths_equidistant_knot_seq(knots, NULL, 0, (size_t)*num_knot_seq, knots, NULL);
}

// Vector math
void ps_vec2_init_(double *out, double *x, double *y) { ts_vec2_init(out, *x, *y); }
void ps_vec3_init_(double *out, double *x, double *y, double *z) { ts_vec3_init(out, *x,*y,*z); }
void ps_vec4_init_(double *out, double *x, double *y, double *z, double *w) { ts_vec4_init(out,*x,*y,*z,*w); }
void ps_vec_add_(double *x, double *y, int *dim, double *out) { ts_vec_add(x, y, (size_t)*dim, out); }
void ps_vec_sub_(double *x, double *y, int *dim, double *out) { ts_vec_sub(x, y, (size_t)*dim, out); }
double ps_vec_dot_(double *x, double *y, int *dim) { return ts_vec_dot(x, y, (size_t)*dim); }
double ps_vec_mag_(double *x, int *dim) { return ts_vec_mag(x, (size_t)*dim); }
void ps_vec_norm_(double *x, int *dim, double *out) { ts_vec_norm(x, (size_t)*dim, out); }
void ps_vec_cross_(double *x, double *y, double *out) { ts_vec3_cross(x,y,out); }

// 設定變距容差
//int ps_fequals_(double *x, double *y) { return ts_fequals(*x,*y); }
