#include "Filling_Holes.h"
#include "common/array.h"
#include "common/edge.h"
#include "common/triangle.h"
#include <set>
#include <map>

struct WeightSet
{
    double angle = 0;
    double area = 0;

    WeightSet() {}
    WeightSet(double angle, double area) {
        this->angle = angle;
        this->area = area;
    }

    //演算子をオーバーロード
    inline WeightSet operator+(const WeightSet& ws) {
        WeightSet res;
        res.angle = std::max(this->angle, ws.angle);
        res.area = this->area + ws.area;
        return res;
    }

    inline const WeightSet operator+(const WeightSet& ws) const {
        WeightSet res;
        res.angle = std::max(this->angle, ws.angle);
        res.area = this->area + ws.area;
        return res;
    }

    //比較演算子をオーバーロード
    bool operator==(const WeightSet rhs) {
        return this->angle == rhs.angle && this->area == rhs.area;
    }
    const bool operator==(const WeightSet rhs) const {
        return this->angle == rhs.angle && this->area == rhs.area;
    }

    bool operator!=(const WeightSet rhs) {
        return this->angle != rhs.angle || this->area != rhs.area;
    }
    const bool operator!=(const WeightSet rhs) const {
        return this->angle != rhs.angle || this->area != rhs.area;
    }

    bool operator<(const WeightSet rhs) {
        return this->angle < rhs.angle || (this->angle == rhs.angle && this->area < rhs.area);
    }
    const bool operator<(const WeightSet rhs) const {
        return this->angle < rhs.angle || (this->angle == rhs.angle && this->area < rhs.area);
    }

    bool operator<=(const WeightSet rhs) {
        return this->angle <= rhs.angle || (this->angle == rhs.angle && this->area <= rhs.area);
    }
    const bool operator<=(const WeightSet rhs) const {
        return this->angle <= rhs.angle || (this->angle == rhs.angle && this->area <= rhs.area);
    }

    bool operator>(const WeightSet rhs) {
        return this->angle > rhs.angle || (this->angle == rhs.angle && this->area > rhs.area);
    }
    const bool operator>(const WeightSet rhs) const {
        return this->angle > rhs.angle || (this->angle == rhs.angle && this->area > rhs.area);
    }

    bool operator>=(const WeightSet rhs) {
        return this->angle >= rhs.angle || (this->angle == rhs.angle && this->area >= rhs.area);
    }
    const bool operator>=(const WeightSet rhs) const {
        return this->angle >= rhs.angle || (this->angle == rhs.angle && this->area >= rhs.area);
    }
};

//最適な三角形を格納して返す再帰関数
void extractDivideTriangle(const int i, const int k, std::vector<Triangle>& dividedTriangles, const std::vector<std::vector<int>>& lambda, const std::vector<int>& point_idx_at_Hole, const std::vector<Point>& points) {
    Point pointI = points[point_idx_at_Hole[i]];
    Point pointI_1 = points[point_idx_at_Hole[i + 1]];
    Point pointK = points[point_idx_at_Hole[k]];
    if (i + 2 == k) dividedTriangles.push_back(Triangle(pointI, pointI_1, pointK));
    else {
        int o = lambda[i][k];
        Point pointO = points[point_idx_at_Hole[o]];
        if (o != i + 1) extractDivideTriangle(i, o, dividedTriangles, lambda, point_idx_at_Hole, points);
        dividedTriangles.push_back(Triangle(pointI, pointO, pointK));
        if (o != k - 1) extractDivideTriangle(o, k, dividedTriangles, lambda, point_idx_at_Hole, points);
    }
}

//3つ以上の辺を持つ頂点を有する三角形をgenerateTrianglesから削除する
void deleteTriangles(std::vector<Triangle>& generateTriangles, const std::vector<Point>& points) {
    //一旦setに入れて削除し、vectorに戻すことで処理を簡単にする
    std::set<Triangle> generateTri_set;
    for (int i = 0; i < generateTriangles.size(); i++) {
        generateTri_set.insert(generateTriangles[i]);
    }

    //point2tri, boundaryedge構築
    std::vector<std::vector<Triangle>> point2tri(points.size(), std::vector<Triangle>(0));
    std::set<Edge> boundaryEdge;
    for (int i = 0; i < generateTriangles.size(); i++) {
        //三角形から頂点のインデックスを取り出し、それを使って辺を作る
        const int ver_idx1 = generateTriangles[i].vertex1.index;
        const int ver_idx2 = generateTriangles[i].vertex2.index;
        const int ver_idx3 = generateTriangles[i].vertex3.index;
        Edge e1 = Edge(points[ver_idx1], points[ver_idx2]);
        Edge e2 = Edge(points[ver_idx2], points[ver_idx3]);
        Edge e3 = Edge(points[ver_idx3], points[ver_idx1]);

        //point2tri構築
        Triangle tri = generateTriangles[i];
        point2tri[ver_idx1].push_back(tri);
        point2tri[ver_idx2].push_back(tri);
        point2tri[ver_idx3].push_back(tri);

        //boundaryEdge構築
        if (boundaryEdge.find(e1) == boundaryEdge.end()) boundaryEdge.insert(e1);
        else boundaryEdge.erase(e1);
        if (boundaryEdge.find(e2) == boundaryEdge.end()) boundaryEdge.insert(e2);
        else boundaryEdge.erase(e2);
        if (boundaryEdge.find(e3) == boundaryEdge.end()) boundaryEdge.insert(e3);
        else boundaryEdge.erase(e3);
    }

    //point2Edgeを構築
    std::vector<std::vector<Edge>> point2edge(points.size(), std::vector<Edge>(0));
    for (auto e : boundaryEdge) {
        point2edge[e.vertex1.index].push_back(e);
        point2edge[e.vertex2.index].push_back(e);
    }

    for (int p_idx = 0; p_idx < point2edge.size(); p_idx++) {
        if (point2edge[p_idx].size() > 2) {
            //該当する辺から三角形を探索してこれらの三角形を削除する
            for (int j = 0; j < point2tri[p_idx].size(); j++) {
                generateTri_set.erase(point2tri[p_idx][j]);
            }
        }
    }

    //削除処理を行って残った三角形をgenerateTrianglesに入れる
    generateTriangles.clear();
    for (auto tri : generateTri_set) {
        generateTriangles.push_back(tri);
    }
}

void fillByWeightFunction(std::vector<Point>& points, std::vector<Triangle>& generateTriangles)
{
    int previousSized = -1;
    while (previousSized != generateTriangles.size()) {
        previousSized = generateTriangles.size();
        deleteTriangles(generateTriangles, points);
    }
    printf("削除後の三角形: %d\n", generateTriangles.size());

    //穴埋めの操作を更新がなくなるまで繰り返す
    int previousSize = -1;

    while (previousSize != generateTriangles.size()) {
        //水密なメッシュにするために分岐点が3つ以上の点を有する三角形を削除する
        previousSize = generateTriangles.size();

        std::set<Edge> boundaryEdge;
        std::vector<std::vector<Edge>> point2edge(points.size(), std::vector<Edge>(0));
        std::map<Edge, Triangle> edge2tri;

        for (int i = 0; i < generateTriangles.size(); i++) {
            //三角形から頂点のインデックスを取り出し、それを使って辺を作る
            const int ver_idx1 = generateTriangles[i].vertex1.index;
            const int ver_idx2 = generateTriangles[i].vertex2.index;
            const int ver_idx3 = generateTriangles[i].vertex3.index;

            Edge e1 = Edge(points[ver_idx1], points[ver_idx2]);
            Edge e2 = Edge(points[ver_idx2], points[ver_idx3]);
            Edge e3 = Edge(points[ver_idx3], points[ver_idx1]);

            //同時にedge2triを構築
            Triangle tri(points[ver_idx1], points[ver_idx2], points[ver_idx3]);
            edge2tri.insert(std::make_pair(e1, tri)); edge2tri.insert(std::make_pair(e2, tri)); edge2tri.insert(std::make_pair(e3, tri));

            //setに辺を入れる、もしくは削除する
            if (boundaryEdge.find(e1) == boundaryEdge.end()) boundaryEdge.insert(e1);
            else boundaryEdge.erase(e1);
            if (boundaryEdge.find(e2) == boundaryEdge.end()) boundaryEdge.insert(e2);
            else boundaryEdge.erase(e2);
            if (boundaryEdge.find(e3) == boundaryEdge.end()) boundaryEdge.insert(e3);
            else boundaryEdge.erase(e3);
        }
        printf("境界辺の数(穴埋め前): %d\n", boundaryEdge.size());

        for (auto e : boundaryEdge) {
            //point2Edgeを構築
            point2edge[e.vertex1.index].push_back(e);
            point2edge[e.vertex2.index].push_back(e);
        }

        while (boundaryEdge.size() != 0) {
            //処理を行う辺を取得、setから削除
            std::set<Edge>::iterator itr = boundaryEdge.begin();
            Edge e = *itr;
            boundaryEdge.erase(e);

            //トレース時のスタートとゴールとなる点を定義して、穴を構成する点を入れる配列を用意
            const int startPointIdx = e.vertex1.index;
            int goalPointIdx = e.vertex2.index;
            std::vector<int> point_idx_at_Hole;
            point_idx_at_Hole.push_back(startPointIdx);

            //穴を構成する全ての点を見つける
            int nowPointIdx = startPointIdx;
            Edge nowEdge = e;
            int cnt = 0;
            bool reverseTrace = false;
            while (nowPointIdx != goalPointIdx) {
                //点が2つの辺を構成していないときは反対周りにトレースする
                int point_has_edge_num = point2edge[nowPointIdx].size();
                if (point_has_edge_num != 2) {
                    if (reverseTrace) break;
                    int tmp = goalPointIdx;
                    goalPointIdx = nowPointIdx;
                    nowPointIdx = tmp;
                    std::reverse(point_idx_at_Hole.begin(), point_idx_at_Hole.end());
                    point_idx_at_Hole.push_back(nowPointIdx);
                    reverseTrace = true;
                    continue;
                }

                //点によって辺を移動
                if (boundaryEdge.find(point2edge[nowPointIdx][0]) != boundaryEdge.end()) nowEdge = point2edge[nowPointIdx][0];
                else nowEdge = point2edge[nowPointIdx][1];

                //辺によって点を移動
                if (nowEdge.vertex1.index == nowPointIdx) nowPointIdx = nowEdge.vertex2.index;
                else nowPointIdx = nowEdge.vertex1.index;

                //穴を構成する点を配列に入れ、setから辺を削除
                point_idx_at_Hole.push_back(nowPointIdx);
                boundaryEdge.erase(nowEdge);

                cnt++;
                if (cnt > 150) break;
            }

            if (point_idx_at_Hole[0] == point_idx_at_Hole[point_idx_at_Hole.size() - 1]) {
                point_idx_at_Hole.resize(point_idx_at_Hole.size() - 1);
            }

            const int pN = point_idx_at_Hole.size();

            //穴が小さすぎたり探索が途中で終わっていたらcontinueする
            if (nowPointIdx != goalPointIdx) continue;
            if (pN <= 2) continue;
            if (pN == 3) {
                auto searchResult = std::find(generateTriangles.begin(), generateTriangles.end(), Triangle(points[point_idx_at_Hole[0]], points[point_idx_at_Hole[1]], points[point_idx_at_Hole[2]]));
                if (searchResult != generateTriangles.end()) continue;
            }

            //穴を埋めるための変数の準備(最初に全ての値を初期化。これによって対角成分の初期化をスキップ)
            WeightSet ws0;
            std::vector<std::vector<WeightSet>> Weights(pN, std::vector<WeightSet>(pN, ws0));
            std::vector<std::vector<int>> lambda(pN - 1, std::vector<int>(pN, 0));

            //隣接する三角形との最大二面角と三角形の面積を初期値として代入する
            for (int i = 0; i < pN - 3; i++) {
                //注目する三角形の点
                Point s = points[point_idx_at_Hole[i]];
                Point t = points[point_idx_at_Hole[i + 1]];
                Point u = points[point_idx_at_Hole[i + 2]];

                //点から辺、辺から隣接する三角形へと探索する。[i, i+2]の辺は存在しないのでここでinsertするとedge2triでエラーを吐く。
                std::set<Edge> holeEdge;
                holeEdge.insert(Edge(s, t));
                holeEdge.insert(Edge(t, u));

                //隣接する三角形を探索すると同時に飛び出たメッシュを索敵する
                std::vector<Triangle> adjacentTri(0);
                for (auto e : holeEdge) {
                    Triangle tri = edge2tri.at(e);
                    adjacentTri.push_back(tri);
                }

                //隣接三角形との二面角を求めていき、最も大きいものを保持
                Triangle holeTri = Triangle(s, t, u);   //注目している三角形
                Array holeTriNormal = holeTri.calcNormal();
                double maxAngle = -1.0;
                for (int j = 0; j < adjacentTri.size(); j++) {
                    Array adjaTriNormal = adjacentTri[j].calcNormal();
                    double angle = holeTriNormal.calcAngle(adjaTriNormal);
                    if (angle < maxAngle) {
                        maxAngle = angle;
                    }
                }

                //面積を計算し、二面角とともに初期値として格納
                double area = holeTri.calcArea();
                WeightSet ws(maxAngle, area);
                Weights[i][i + 2] = ws;

                //論文には書いていないが、こうするしかないのでlambdaも初期化
                lambda[i][i + 1] = 0;   //ここ合ってるか怪しい？
                lambda[i][i + 2] = i + 1;
            }

            //動的計画法でWeightsを埋めていく。最小のmをlambdaに入れる
            for (int j = 3; j < pN; j++) {
                for (int i = 0; i < pN - j; i++) {
                    //使用する点や添え字
                    int k = i + j;
                    Point pointI = points[point_idx_at_Hole[i]];
                    Point pointK = points[point_idx_at_Hole[k]];

                    //重みが最小になるようなmを求める
                    int minM = i + 1;
                    WeightSet minWS = WeightSet(100, 100);
                    for (int m = i + 1; m < k; m++) {
                        Point pointM = points[point_idx_at_Hole[m]];
                        Triangle usingMTri = Triangle(pointI, pointM, pointK);  //点Mを使った三角形
                        Array usingM_normal = usingMTri.calcNormal();

                        //二面角を求める。4点の穴の場合は論文通りだと三角形が作れない場合が存在するため、三角形が作れるかどうかで分岐させる
                        double dihedral = 0;

                        if (i + 1 == m) {
                            Triangle triM2K = Triangle(pointM, points[point_idx_at_Hole[lambda[m][k]]], pointK);
                            Array M2K_normal = triM2K.calcNormal();
                            dihedral = usingM_normal.calcAngle(M2K_normal);
                        }
                        else if (m + 1 == k) {
                            Triangle triI2M = Triangle(pointI, points[point_idx_at_Hole[lambda[i][m]]], pointM);
                            Array I2M_noraml = triI2M.calcNormal();
                            dihedral = usingM_normal.calcAngle(I2M_noraml);
                        }
                        else {
                            Triangle triI2M = Triangle(pointI, points[point_idx_at_Hole[lambda[i][m]]], pointM);
                            Triangle triM2K = Triangle(pointM, points[point_idx_at_Hole[lambda[m][k]]], pointK);
                            Array I2M_noraml = triI2M.calcNormal();
                            Array M2K_normal = triM2K.calcNormal();
                            double dihedral_mk = usingM_normal.calcAngle(M2K_normal);
                            double dihedral_im = usingM_normal.calcAngle(I2M_noraml);
                            dihedral = std::max(dihedral_im, dihedral_mk);
                        }

                        //最後はi,kが隣接辺になるのでさらに比較する
                        if (i == 0 && k == pN - 1) {
                            Triangle triI2K = edge2tri.at(Edge(points[point_idx_at_Hole[i]], points[point_idx_at_Hole[k]]));
                            Array I2K_normal = triI2K.calcNormal();
                            double dihedral_ik = usingM_normal.calcAngle(I2K_normal);
                            dihedral = std::max(dihedral, dihedral_ik);
                        }

                        //面積を求める
                        double area = usingMTri.calcArea();
                        WeightSet ws = WeightSet(Weights[i][m] + Weights[m][k] + WeightSet(dihedral, area));

                        //最小値を更新
                        if (ws < minWS) {
                            minWS = ws;
                            minM = m;
                        }
                    }
                    //配列を更新
                    Weights[i][k] = minWS;
                    lambda[i][k] = minM;
                }
            }
            //最適な三角形分割を抽出する
            std::vector<Triangle> divededTri;
            extractDivideTriangle(0, pN - 1, divededTri, lambda, point_idx_at_Hole, points);

            //抽出した三角形をgenerateTrianglesに追加する
            for (int i = 0; i < divededTri.size(); i++) {
                if (divededTri[i].vertex1 == divededTri[i].vertex2 || divededTri[i].vertex2 == divededTri[i].vertex3 || divededTri[i].vertex3 == divededTri[i].vertex1) {
                    continue;
                }
                generateTriangles.push_back(Triangle(divededTri[i].vertex1, divededTri[i].vertex2, divededTri[i].vertex3));
            }
        }
    }
}