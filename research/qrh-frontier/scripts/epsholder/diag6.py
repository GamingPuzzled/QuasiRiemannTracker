from cond import theta_float, sixth, geom, feasible
for e in [1.1801, 1.18, 1.175, 1.17, 1.16, 1.15, 1.10, 1.0]:
    fl = theta_float(sixth(e))
    l, h = fl['l'], fl['h']; x = 1+l-h; y = 1-l-x
    print(f"e6={e}: theta={fl['theta']:.13f} l={l:.10f} h={h:.10f} x={x:.6f} y={y:.6f} y-x={y-x:.2e} "
          f"l/h-4/21={l/h-4/21:.2e} 4l-y={4*l-y:.2e} 2x-y-l={2*x-y-l:.2e} delta*={fl['delta']:.5f} act={fl['active']} gap={fl['gap']:.1e}")
