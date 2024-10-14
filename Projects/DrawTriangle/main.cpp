#include <tgaimage.h>
#include <model.h>
#include <glm.hpp>

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

const TGAColor white = { 255, 255, 255, 255 };
const TGAColor red = { 0, 0, 255, 255 };
const TGAColor blue = { 255, 0, 0, 255 };

const std::string modelPath = "Res/african_head.obj";
const std::string diffuseTexture = "Res/african_head_diffuse.tga";
const int width = 200;
const int height = 200;

void line(glm::vec2 vert0, glm::vec2 vert1, TGAImage& image, TGAColor color)
{
	int x0 = (int)vert0.x;
	int y0 = (int)vert0.y;
	int x1 = (int)vert1.x;
	int y1 = (int)vert1.y;
	// print x0, y0, x1, y1
	std::cout << "(" << x0 << ", " << y0 << ") -> (" << x1 << ", " << y1 << ")\n";

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

//void Triangle(glm::vec2 vert0, glm::vec2 vert1, glm::vec2 vert2, TGAImage& image, TGAColor color)
//{
//    line(vert0, vert1, image, color);
//	line(vert1, vert2, image, color);
//	line(vert2, vert0, image, color);
//}

void barycentric(glm::vec3* vertices, glm::vec2 P, float& u, float& v, float& w)
{
    glm::vec2 vert0 = glm::vec2(vertices[0].x, vertices[0].y);
    glm::vec2 vert1 = glm::vec2(vertices[1].x, vertices[1].y);
    glm::vec2 vert2 = glm::vec2(vertices[2].x, vertices[2].y);
    glm::vec2 ab = vert1 - vert0; // AB
    glm::vec2 ac = vert2 - vert0; // AC
    glm::vec2 ap = P - vert0; // AP

    float d00 = glm::dot(ab, ab); // |AB|^2
    float d01 = glm::dot(ab, ac); 
    float d11 = glm::dot(ac, ac); // |AC|^2
    float d20 = glm::dot(ap, ab);
    float d21 = glm::dot(ap, ac);
    float denom = 1. / (d00 * d11 - d01 * d01); 
    v = (d11 * d20 - d01 * d21) * denom;
    w = (d00 * d21 - d01 * d20) * denom;
    u = 1. - v - w;
}

void Triangle(glm::vec3* vertices, glm::vec2* uvs, float* zBuffer, TGAImage& image, TGAImage& diffusemap)
{
    // Calculate Bounding Box
	int minX = std::min(std::min(vertices[0].x, vertices[1].x), vertices[2].x);
    int maxX = std::max(std::max(vertices[0].x, vertices[1].x), vertices[2].x);
    int minY = std::min(std::min(vertices[0].y, vertices[1].y), vertices[2].y);
    int maxY = std::max(std::max(vertices[0].y, vertices[1].y), vertices[2].y);

    int width = image.width()  - 1;
	int height = image.height() - 1;
    minX = std::max(0, minX);
    minY = std::max(0, minY);
    maxX = std::min(width, maxX);
	maxY = std::min(height, maxY);

    glm::vec3 p;
	// Iterate over each pixel in bounding box
	for (p.x = minX; p.x <= maxX; p.x++)
	{
		for (p.y = minY; p.y <= maxY; p.y++)
		{
            float u, v, w;
            barycentric(vertices, p, u, v, w);
            glm::vec3 bc = glm::vec3(u, v, w);
            p.z = 0.;
            if (u < 0 || v < 0 || w < 0)
                continue;
			// Interpolate z
            for (int i = 0; i < 3; i++) 
            {
                p.z += vertices[i].z * bc[i];
            }

			// Check if pixel is inside triangle
			if (zBuffer[int(p.x + p.y * width)] < p.z)
			{
				zBuffer[int(p.x + p.y * width)] = p.z;

                // Interpolate uv
                glm::vec2 uv = glm::vec2(0, 0);
                for (int i = 0; i < 3; i++) 
                {
                    uv.x += uvs[i].x * bc[i];
                    uv.y += uvs[i].y * bc[i];
                }

                // Get diffuse color from diffuse texture
                TGAColor color = diffusemap.get((int)(uv.x * diffusemap.width() - 0.5f), (int)(uv.y * diffusemap.height() - 0.5f));
                
                // output depth value
				/*float depth = (int)((p.z * 0.5f + 0.5f) * 255);
                color.bgra[0] = depth;
                color.bgra[1] = depth;
                color.bgra[2] = depth;
                color.bgra[3] = 255;*/
                
				// Draw pixel on image
				image.set(p.x, p.y, color);
			}
		}
	}
}

int main(int argc, char** argv) {
    TGAImage image(width, height, TGAImage::RGB);

    //glm::vec2 triangle0[3] = { glm::vec2(10, 10), glm::vec2(100, 30), glm::vec2(190, 160) };

	//Triangle(triangle0, image, red);
  /*  glm::vec2 t0[3] = {glm::vec2(10, 70),   glm::vec2(50, 160),  glm::vec2(70, 80)}; 
    glm::vec2 t1[3] = {glm::vec2(180, 50),  glm::vec2(150, 1),   glm::vec2(70, 180)}; 
    glm::vec2 t2[3] = {glm::vec2(180, 150), glm::vec2(120, 160), glm::vec2(130, 180)}; */

  /*  Triangle(t0[0], t0[1], t0[2], image, red);
	Triangle(t1[0], t1[1], t1[2], image, white);
	Triangle(t2[0], t2[1], t2[2], image, blue);*/

    Model model(modelPath);
    glm::vec3 lightDir = glm::vec3(0, 0, -1);
    float* zBuffer = new float[width * height];

    // load diffuse texture
    //stbi_set_flip_vertically_on_load(true);
    int texWidth, texHeight, texChannels;
    stbi_uc* pixels = stbi_load(diffuseTexture.c_str(), &texWidth, &texHeight, &texChannels, STBI_rgb_alpha);
    TGAImage diffusemap(texWidth, texHeight, TGAImage::RGBA);
    if (pixels)
    {
        // create diffuse texture
        for (int y = 0; y < texHeight; y++)
        {
            for (int x = 0; x < texWidth; x++)
            {
                TGAColor color;
                color.bgra[2] = pixels[(y * texWidth + x) * 4 + 0];
                color.bgra[1] = pixels[(y * texWidth + x) * 4 + 1];
                color.bgra[0] = pixels[(y * texWidth + x) * 4 + 2];
                color.bgra[3] = pixels[(y * texWidth + x) * 4 + 3];
                diffusemap.set(x, y, color);
            }
        }

        stbi_image_free(pixels);
    }
    else
    {
        std::cout << "Failed to load texture: " << diffuseTexture << std::endl;
    }
    
    for (int i = 0; i < model.nfaces(); i++)
    {
        // verteces
        glm::vec3 vertices[3];
        // trasform vertices to screen space
        glm::vec3 vertices2[3];
        // colors
        // TGAColor colors[3];
        
        // uvs
        glm::vec2 uvs[3];
        
        for (int j = 0; j < 3; j++)
        {
            vec3 vertex = model.vert(i, j);            
			vertices2[j] = glm::vec3(vertex.x, vertex.y, vertex.z);
            vertices[j] = glm::vec3((vertex.x + 1.0) * width / 2, (vertex.y + 1.0) * height / 2, vertex.z);

            // calculate color
            // colors[j] = diffusemap.get((int)(model.uv(i, j).x * texWidth - 0.5f), (int)(model.uv(i, j).y * texHeight - 0.5f));

            // calculate uv
            uvs[j] = glm::vec2(model.uv(i, j).x, model.uv(i, j).y);
        }

        // calculate normal
        glm::vec3 normal = glm::cross(vertices2[2] - vertices2[0], vertices2[1] - vertices2[0]);
        normal = glm::normalize(normal);
		float intensity = glm::dot(lightDir, normal);
        if (intensity > 0)
        {
            // TGAColor color = { intensity * 255, intensity * 255, intensity * 255, 255 };
            Triangle(vertices, uvs, zBuffer, image, diffusemap);
        }
       
    }
            
    //image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
    image.write_tga_file("output.tga");
    return 0;
}