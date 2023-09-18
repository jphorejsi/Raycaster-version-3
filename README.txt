directory for glm.hpp needs to installed to this folder.

use "make" to prompt program and ./target "text-file" to run the program, a new .ppm file will be created that will resemble the conditions of the text file.

text file rules:
imsize: resolution of the outputted .ppm file.
eye: location of the camera/eye we are seeing from in 3d space.
viewdir: the viewing direction, the direction the camera is facing.
hfov: horizontal fov.
updir: the upwards direction of the eye/camera, can be used to rotate what is being viewed.
bkgcolor: the color to be rendered if a ray does not intersect an object, the background's color.
mtlcolor: properties for the coloring and shading of the objects referenced to after. Includes object's diffuse color, ambient coefficient, diffuse coefficient, specular coefficent and specular exponent.
sphere: the location and radius of a sphere in 3d space.
light: the location of the light in 3d space, the type of light it is (1 for point, 2 for directional), and the color of the light.
v: list of vertex positions
f: list of triangle definitions
vn: list of vertex normal vectors

enables the user to specify, for each smooth shaded triangle, indices into the array of vertex normal directions as well as into the array of vertex positions. Please use this syntax to define smooth-shaded triangles:
  f   v1//vn1    v2//vn2    v3//vn3

In this example, v1, v2, v3 are indices into the array of vertex locations and vn1, vn2, vn3 are indices into the array of normal directions.