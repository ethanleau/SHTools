
#include "SH.h"

#include <iostream>
#include <string>
#include <vector>
#include <cassert>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

class Image
{
public:
	Image(int w, int h, int d, float* dat) : width(w), depth(d), data(dat) { assert(w == h); }

	float& at(int i, int j, int k)
	{
		return data[(i * width + j) * depth + k];
	}

	float at(int i, int j, int k) const
	{
		return data[(i * width + j) * depth + k];
	}

	int width, depth;

private:
	float* data;
};


float sinc(float x) 
{               /* Supporting sinc function */
	if (fabs(x) < 1.0e-4) return 1.0;
	else return(sin(x) / x);
}


// Refer to the paper: 2002 An Efficient Representation for Irradiance Environment Maps
bool sphereMapping(int i, int j, int width, float& x, float& y, float& z, float& domega)
{
	float u, v, r, theta, phi;

	v = (width / 2.0 - i) / (width / 2.0);		// v ranges from -1 to 1
	u = (j - width / 2.0) / (width / 2.0);		// u ranges from -1 to 1
	r = sqrt(u * u + v * v);

	// Consider only circle with r < 1
	if (r > 1.0)
		return false;

	theta = PI * r;
	phi = atan2(v, u);

	// Cartesian components
	x = sin(theta) * cos(phi);
	y = sin(theta) * sin(phi);
	z = cos(theta);

	// Computation of the solid angle.  This follows from some elementary calculus converting sin(theta) d theta d phi into
    // coordinates in terms of r.  This calculation should be redone if the form of the input changes
	domega = (2 * PI / width) * (2 * PI / width) * sinc(theta);

	return true;
}


SHColor project(FILTER_TYPE filterType, const Image& img)
{
	SHColor shColor;

	for (int i = 0; i < img.width; i++)
	{
		for (int j = 0; j < img.width; j++)
		{
			// We now find the cartesian components for the point (i,j)
			float x, y, z, domega;
			if (sphereMapping(i, j, img.width, x, y, z, domega))
			{
				float color[3] = { img.at(i, j, 0), img.at(i, j, 1), img.at(i,j,2) };

				// Update Integration
				SHCoefficients shDirection;
				projectOnSH(x, y, z, shDirection);

				switch (filterType)
				{
				case FILTER_GAUSSIAN:
					gauss(SH_ORDER, shDirection);
					break;
				case FILTER_HANNING:
					hanning(SH_ORDER, shDirection);
					break;
				case FILTER_LANCZOS:
					lanczos(SH_ORDER, shDirection);
					break;
				default:
					break;
				}

				for (int c = 0; c < 3; ++c)
				{
					for (int k = 0; k < SH_COUNT; ++k)
					{
						shColor[c][k] += shDirection[k] * color[c] * domega;
					}
				}
			}
		}
	}
	return shColor;
}

void reconstruct(Image& img, const SHColor& shColor)
{
	for (int i = 0; i < img.width; i++)
	{
		for (int j = 0; j < img.width; j++)
		{
			// We now find the cartesian components for the point (i,j)
			float x, y, z, domega;
			if (sphereMapping(i, j, img.width, x, y, z, domega))
			{
				float color[3];
				reconstructFromSH(x, y, z, shColor, color);

				for (int c = 0; c < 3; ++c)
					img.at(i, j, c) = color[c];
			}
		}
	}
}

void printUsageAndExit(const char* argv0)
{
	std::cerr << "Usage  : " << argv0 << " [options] <input>\n";
	std::cerr << "<input>              .hdr probe image\n";
	std::cerr << "Options: --help      Print this usage message\n";
	std::cerr << "         -h          Hanning  filter\n";
	std::cerr << "         -l          Lanczos  filter\n";
	std::cerr << "         -g          Gaussian filter\n";
	exit(0);
}

int main(int argc, char** argv)
{
	FILTER_TYPE filterType = FILTER_DISABLE;
	std::string infile;

	for (int i = 1; i < argc; ++i)
	{
		const std::string arg = argv[i];
		if (arg == "--help")
		{
			printUsageAndExit(argv[0]);
		}
		else if (arg == "-l")
		{
			filterType = FILTER_LANCZOS;
		}
		else if (arg == "-h")
		{
			filterType = FILTER_HANNING;
		}
		else if (arg == "-g")
		{
			filterType = FILTER_GAUSSIAN;
		}
		else
		{
			infile = arg;
		}
	}

	if (infile.empty())
	{
		std::cerr << "need input file argument" << std::endl;
		printUsageAndExit(argv[0]);
	}

	std::string outfile = infile.substr(0, infile.size() - 4) + "-" + std::to_string(SH_ORDER);
	switch (filterType)
	{
	case FILTER_DISABLE:
		break;
	case FILTER_GAUSSIAN:
		outfile += "-g";
		break;
	case FILTER_HANNING:
		outfile += "-h";
		break;
	case FILTER_LANCZOS:
		outfile += "l";
		break;
	default:
		break;
	}
	outfile += ".hdr";

	int w, h, d;

	// Read .HDR image
	float* data = stbi_loadf(infile.c_str(), &w, &h, &d, 0);
	if (data == nullptr)
	{
		std::cerr << "input file argument is incorrect" << std::endl;
		printUsageAndExit(argv[0]);
	}

	Image img(w, h, d, data);

	// Project on SH
	SHColor shColor = project(filterType, img);

	// Log
	for (int i = 0; i < SH_COUNT; ++i)
	{
		printf_s("m_params.env_radiance[%d] = make_float3(%9.6f, %9.6f, %9.6f);\n", i, shColor[0][i], shColor[1][i], shColor[2][i]);
	}

	// Reconstruct
	float* newData = new float[w * h * d];
	Image newImg(w, h, d, newData);
	reconstruct(newImg, shColor);
	stbi_write_hdr(outfile.c_str(), w, h, d, newData);

	stbi_image_free(data);
	delete newData;
	return 0;
}
