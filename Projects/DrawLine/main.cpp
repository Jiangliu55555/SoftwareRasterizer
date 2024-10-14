#include <tgaimage.h>
#include <model.h>

const TGAColor white = { 255, 255, 255, 255 };
const TGAColor red = { 0, 0, 255, 255 };
const std::string modelPath = "Res/african_head.obj";
const int width = 1024;
const int height = 1024;

void line(int x0, int y0, int x1, int y1, TGAImage& image, TGAColor color)
{
    bool steep = false;
    // low steep axis as interpolation value
    if (std::abs(x1 - x0) < std::abs(y1 - y0))
    {
		std::swap(x0, y0);
		std::swap(x1, y1);
		steep = true;
    }

    // left to right
    if (x0 > x1)
    {
		std::swap(x0, x1);
		std::swap(y0, y1);
    }

    int dx = x1 - x0;
    int dy = y1 - y0;
    int derror = std::abs(dy) * 2;
    int error = 0;
    int y = y0;

    for (int x = x0; x <= x1; x++)
    {
        if (steep)
        {
            image.set(y, x, color);
        }
        else {
            image.set(x, y, color);
        }

        error += derror;
        if (error > dx) {
            y += (y1 > y0 ? 1 : -1);
            error -= dx * 2;
        }
    }
}

int main(int argc, char** argv) {
    TGAImage image(width, height, TGAImage::RGB);

    /* for (size_t i = 0; i < 1000000; i++)
     {
         line(0, 0, 50, 50, image, red);
     }*/

     /* line(13, 20, 80, 40, image, white);
      line(20, 13, 40, 80, image, red);
      line(80, 40, 13, 20, image, red);*/
    Model model(modelPath);
    for (int i = 0; i < model.nfaces(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            vec3 v0 = model.vert(i, j);
            vec3 v1 = model.vert(i, (j + 1) % 3);
			// print v0, v1
			//std::cout << "(" << v0.x << ", " << v0.y << ", " << v0.z << ") -> (" << v1.x << ", " << v1.y << ", " << v1.z << ")\n";
            // localPos convert to imagePixel Pos: [-1, 1] -> [0, width/height]
	        int x0 = (v0.x + 1.) * width / 2.;
            int y0 = (v0.y + 1.) * height / 2.;
            int x1 = (v1.x + 1.) * width / 2.;
            int y1 = (v1.y + 1.) * height / 2.;
            line(x0, y0, x1, y1, image, white);
        }
    }
            
    //image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("output.tga");
    return 0;
}