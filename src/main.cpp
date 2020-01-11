#include <iostream>
#include <cmath>
#include <utility>
#include <vector>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <string>
#include "sarge.h"
using namespace std;

typedef std::pair<double, double> Point;

double PerpendicularDistance(const Point &pt, const Point &lineStart, const Point &lineEnd)
{
	double dx = lineEnd.first - lineStart.first;
	double dy = lineEnd.second - lineStart.second;

	//Normalise
	double mag = pow(pow(dx,2.0)+pow(dy,2.0),0.5);
	if(mag > 0.0)
	{
		dx /= mag; dy /= mag;
	}

	double pvx = pt.first - lineStart.first;
	double pvy = pt.second - lineStart.second;

	//Get dot product (project pv onto normalized direction)
	double pvdot = dx * pvx + dy * pvy;

	//Scale line direction vector
	double dsx = pvdot * dx;
	double dsy = pvdot * dy;

	//Subtract this from pv
	double ax = pvx - dsx;
	double ay = pvy - dsy;

	return pow(pow(ax,2.0)+pow(ay,2.0),0.5);
}

void RamerDouglasPeucker(const vector<Point> &pointList, double epsilon, vector<Point> &out)
{
	if(pointList.size()<2)
		throw invalid_argument("Not enough points to simplify");

	// Find the point with the maximum distance from line between start and end
	double dmax = 0.0;
	size_t index = 0;
	size_t end = pointList.size()-1;
	for(size_t i = 1; i < end; i++)
	{
		double d = PerpendicularDistance(pointList[i], pointList[0], pointList[end]);
		if (d > dmax)
		{
			index = i;
			dmax = d;
		}
	}

	// If max distance is greater than epsilon, recursively simplify
	if(dmax > epsilon)
	{
		// Recursive call
		vector<Point> recResults1;
		vector<Point> recResults2;
		vector<Point> firstLine(pointList.begin(), pointList.begin()+index+1);
		vector<Point> lastLine(pointList.begin()+index, pointList.end());
		RamerDouglasPeucker(firstLine, epsilon, recResults1);
		RamerDouglasPeucker(lastLine, epsilon, recResults2);

		// Build the result list
		out.assign(recResults1.begin(), recResults1.end()-1);
		out.insert(out.end(), recResults2.begin(), recResults2.end());
		if(out.size()<2)
			throw runtime_error("Problem assembling output");
	}
	else
	{
		//Just return start and end points
		out.clear();
		out.push_back(pointList[0]);
		out.push_back(pointList[end]);
	}
}


/**
 * Reads csv file into table, exported as a vector of vector of doubles.
 * @param inputFileName input file name (full path).
 * @return data as vector of vector of doubles.
 */
vector< vector<double> > parse2DCsvFile(string inputFileName)
{

    vector<vector<double> > data;
    std::ifstream inputFile(inputFileName.c_str());
    int l = 0;

    while (inputFile) {
        l++;
        string s;
        if (!getline(inputFile, s)) break;
        if (s[0] != '#') {
            istringstream ss(s);
            vector<double> record;

            while (ss) {
                string line;
                if (!getline(ss, line, ','))
                    break;
                try {
                    record.push_back(stof(line));
                }
                catch (const std::invalid_argument e) {
                    cout << "NaN found in file " << inputFileName << " line " << l
                         << endl;
                    e.what();
                }
            }

            data.push_back(record);
        }
    }

    if (!inputFile.eof()) {
        cerr << "Could not read file " << inputFileName << "\n";
    }

    return data;
}


int main(int argc, char** argv)
{
	vector<Point> pointList;
	vector<Point> pointListOut;
	std::string inputfile = "points.csv";
    std::string outputfile = "points_reduced.csv";
    std::string epsilon_str = "0.1";

	Sarge sarge;

    sarge.setArgument("h", "help", "Get help.", false);
    sarge.setArgument("i", "input", "The input file", true);
    sarge.setArgument("o", "output", "The output file", true);
    sarge.setArgument("e", "epsilon", "Distance dimension epsilon", true);
    sarge.setDescription("Reduce number of XY coordinates by using Ramer-Douglas-Peucker algorithm.");
    sarge.setUsage("simplify <options>");

    if (!sarge.parseArguments(argc, argv)) {
        std::cerr << "Couldn't parse arguments..." << std::endl;
        return 1;
    }

    //std::cout << "Number of flags found: " << sarge.flagCount() << std::endl;

    if (sarge.exists("help")) {
        sarge.printHelp();
        return 0;
    }


    if (!sarge.getFlag("input", inputfile))
        std::cout << "No input file defined. Using " << inputfile << std::endl;

    if (!sarge.getFlag("output", outputfile))
        std::cout << "No output file defined. Using " << outputfile << std::endl;

    if (!sarge.getFlag("epsilon", epsilon_str))
        std::cout << "No epsilon defined. Using " << epsilon_str << std::endl;

    double epsilon = std::stod(epsilon_str, nullptr);



	vector<vector<double>> data = parse2DCsvFile(inputfile);

	//cout << "Size of data: " << data.size() << endl;

    for (auto l : data) {
        //cout << "Size of l: " << l.size() << endl;

        if (l.size() == 2)
        {
            int i = 0;
            double point_x;
            double point_y;
            for (auto x : l)
            {
                if (i == 0) point_x = x;
                if (i == 1) point_y = x;
                i++;
            }
            pointList.push_back(Point(point_x, point_y));
        }
    }

	RamerDouglasPeucker(pointList, epsilon, pointListOut);

	// Print to screen
	cout << "Result" << endl;
	for(size_t i=0;i< pointListOut.size();i++)
	{
		cout << pointListOut[i].first << "," << pointListOut[i].second << endl;
	}

	// Write to file
    char delimiter = ',';  // by default: ','
    ofstream fout(outputfile);
    for(size_t i=0;i< pointListOut.size();i++)
    {
        fout << pointListOut[i].first << delimiter << pointListOut[i].second << endl;
    }
    fout.close();

	return 0;
}
