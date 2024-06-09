#include "Filling_Holes.h"
#include "common/io.h"

int main()
{
	//ファイルからメッシュを読み取る
	const char input[] = "D:Mesh/MeshWithHoles.ply";
	int pN;
	std::vector<Point> points;
	std::vector<Triangle> triangles;
	read_ply(input, points, triangles, pN);

	//穴埋め
	fillByWeightFunction(points, triangles);

	//メッシュを出力
	const char output[] = "D:Mesh/MeshWithoutHoles.ply";
	write_ply(output, points, triangles);
}