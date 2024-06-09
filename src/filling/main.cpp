#include "Filling_Holes.h"
#include "common/io.h"

int main()
{
	//�t�@�C�����烁�b�V����ǂݎ��
	const char input[] = "D:Mesh/MeshWithHoles.ply";
	int pN;
	std::vector<Point> points;
	std::vector<Triangle> triangles;
	read_ply(input, points, triangles, pN);

	//������
	fillByWeightFunction(points, triangles);

	//���b�V�����o��
	const char output[] = "D:Mesh/MeshWithoutHoles.ply";
	write_ply(output, points, triangles);
}