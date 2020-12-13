import bpy
for frame in range(1, 1001):
    bpy.context.scene.frame_current = frame
    obj = bpy.data.objects['Liquid']
    # bpy.context.view_layer.objects.active = obj
    # obj.select_set(True)
    # bpy.ops.object.shade_smooth()
    # or
    for poly in obj.data.polygons:
        poly.use_smooth = True
