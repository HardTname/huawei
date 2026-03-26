#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

struct Vector2D;
struct Polygon;

using namespace std;
const double EPS = 1e-6;
// 精度控制函数
inline int Sign(double x){
    if(x>EPS)return 1;
    if(x<EPS)return -1;
    return 0;
}
// 实现二维向量的基本操作
struct Vector2D {
    double x, y;

    Vector2D(double x = 0, double y = 0) : x(x), y(y) {}
    Vector2D operator-(const Vector2D& other) const { return Vector2D(x - other.x, y - other.y); }
    Vector2D operator+(const Vector2D& other) const { return Vector2D(x + other.x, y + other.y); }
    Vector2D operator*(double scalar) const { return Vector2D(x * scalar, y * scalar); }//向量乘以实数
    double Dot(const Vector2D& other) const { return x * other.x + y * other.y; }//向量点积
    double Cross(const Vector2D& other) const {return x*other.y-y*other.x;}
    double Length() const { return std::sqrt(x * x + y * y); }//返回模长
    Vector2D Normalize() const//归一化，求方向上的单位向量
    {
        double len = Length();
        if (len == 0) {
            return *this;
        }
        return Vector2D(x / len, y / len);
    }
    Vector2D Perp() const { return Vector2D(-y, x); }//当前向量逆时针旋转90°
};
// 计算三个点构成的两个向量的叉积 (AB x AC)，用于判断点C相对于向量AB的位置关系
inline double CrossProduct(const Vector2D& a, const Vector2D& b, const Vector2D& c) {
    return (b - a).Cross(c - a);
}
// 点到线段最短距离
struct DistanceInfo{
    double distance;
    Vector2D pushVector;// 点到线段最短距离的移动向量
};
// 点p到线段ab的最短距离和对应的移动向量
DistanceInfo PointToSegmentDistance(const Vector2D& p,const Vector2D& a,const Vector2D& b){
    Vector2D ab = b - a;
    Vector2D ap = p - a;
    double ab_len_sq = ab.Dot(ab);// 求线段 ab 长度的平方
    // 线段两个端点重合
    if (Sign(ab_len_sq) == 0) {
        Vector2D p_to_a=p-a;
        return {p_to_a.Length(), p_to_a};
    }
    // 计算投影点在线段上的比例参数 t = (AP · AB) / |AB|^2 ???这里是干嘛的
    double t = ap.Dot(ab) / ab_len_sq;
    t = std::max(0.0,std::min(1.0,t)); // 将t限制到[0,1]范围内，确保投影点在线段上
    Vector2D closestPoint = a +ab *t; // 计算p在线段上的投影点
    Vector2D pushVector = closestPoint - p; // 计算点p到线段ab的最短移动向量
    return {pushVector.Length(),pushVector};
}
struct Polygon {
    std::vector<Vector2D> vertices;

    Polygon() = default;
    Polygon(std::initializer_list<Vector2D> vts) : vertices(vts) {}
    // 这里有问题？
    Vector2D GetCenter() const
    {
        Vector2D center(0, 0);
        if (vertices.empty()) {
            return center;
        }
        for (const auto& v : vertices) {
            center = center + v;
        }
        return center * (1.0 / vertices.size());
    }
    void MoveByVec(const Vector2D& vec) {
        for (auto& v : vertices) {
            v = v + vec;
        }
    }

    // 判断点 p 与多边形的关系，返回 1 表示在内部，0 表示在边界上，-1 表示在外部
    int Contains(const Vector2D& p) const{
        bool isInside = true;
        int n = vertices.size();

        for(int i = 0,j = n-1; i < n; j = i++){
            const Vector2D &v1 = vertices[i], &v2 = vertices[j];
            
            // 判断点 p 是否这条边上
            if(Sign(CrossProduct(v1,v2,p))==0
               && Sign((p-v1).Dot(p-v2))<=0
               && Sign((p-v1).Dot(p-v2))>=0){
                return -1;
            }

            // 射线法：向 x 轴正方向发射射线，计算交点个数
            if(Sign(v1.y-p.y)!=Sign(v2.y-p.y)){
                // 利用相似三角形计算交点的 x 坐标
                double intersectX = (v2.x - v1.x) * (p.y - v1.y) / (v2.y - v1.y) + v1.x;
                // 如果交点在射线右侧，计数
                if(p.x < intersectX){
                    isInside = !isInside;
                }
            }
        }
        return isInside ? 1 : 0;
    }
};
// upd : 3.26
// 传入多边形，得到关于原点对称的多边形
Polygon GetNativePolygon(const Polygon& poly){
    Polygon negPoly;
    negPoly.vertices.reserve(poly.vertices.size());
    for(const auto& v: poly.vertices){
        //直接将 x 和 y 取反，得到关于原点对称的多边形
        negPoly.vertices.emplace_back(-v.x,-v.y);
    }
    return negPoly;
}

// 检查点 p 是否在三角形 ABC 内部，不包含边界
bool IsPointInTiangle(const Vector2D& p,const Vector2D& a,const Vector2D& b,const Vector2D& c){
    double cp1=CrossProduct(a,b,p);
    double cp2=CrossProduct(b,c,p);
    double cp3=CrossProduct(c,a,p);
    // 如果三个叉积的均大于 0 ，则点在三角形内部（题目保证逆时针输入）
    return Sign(cp1)>0 && Sign(cp2)>0 && Sign(cp3)>0;
}

//将传入的简单多边形分解为多个三角形 （凸多边形）（O(n^3)、性能过差，需要优化）
std::vector<Polygon>PartPolygonToTriangle(const Polygon& poly){
    std::vector<Polygon> triangles;
    std::vector<Vector2D> V = poly.vertices;

    //如果已经是三角形，则直接返回
    if(V.size()<=3){
        triangles.emplace_back(poly);
        return triangles;
    }

    // 不断切除，直至只剩三个顶点(耳切法）
    while(V.size()>3){
        int n = V.size();
        bool earFound = false;

        for(int i=0;i<n;i++){
            int prev = (i-1+n)%n;
            int next=(i+1)%n;

            // 判断角 prev-i-next 是否为凸角（叉积prev_i x next_i > 0）
            if(Sign(CrossProduct(V[prev],V[i],V[next]))>0){
                bool isEar = true;
                // 检查多边形内其他点是否在这个三角形内部（凹角对应的顶点可能会在里面）
                for (int j = 0;j<n;j++){
                    if(j==prev || j==i ||j==next){
                        continue;
                    }
                    if(IsPointInTiangle(V[j],V[prev],V[i],V[next])){
                        isEar=false;
                        break;
                    }
                }
                // 如果是耳朵，就切除掉
                if(isEar){
                    triangles.emplace_back(Polygon({V[prev],V[i],V[next]}));
                    V.erase(V.begin()+i);
                    earFound=true;
                    break;
                }
            } 
        }

        if(!earFound)break;

    }

    if(V.size()==3){
        triangles.emplace_back(Polygon({V[0],V[1],V[2]}));
    }
    return triangles;
}

// 重新排列顶点
void Reorder(std::vector<Vector2D>& points){
    int minidx=0;
    int n=points.size();
    for(int i = 1;i<n;i++){
        if(points[i].y<points[minidx].y||(Sign(points[i].y-points[minidx].y)==0 && points[i].x<points[minidx].x)){
            minidx=i;
        }
    }
    std::rotate(points.begin(),points.begin()+minidx,points.end());
}

//计算凸多边形的闵可夫斯基和(有什么用?)
Polygon MinkowskiSum(Polygon A,Polygon B){
    Reorder(A.vertices);
    Reorder(B.vertices);
    
    // 补齐末尾顶点
    A.vertices.push_back(A.vertices[0]);
    B.vertices.push_back(B.vertices[0]);

    Polygon result;
    int i=0,j=0;
    int n = A.vertices.size()-1;
    int m = B.vertices.size()-1;

    while(i<n||j<m){
        result.vertices.push_back(A.vertices[i]+B.vertices[i]);

        double cross=0;
        if(i<n&&j<m){
            Vector2D edgeA = A.vertices[i+1]-A.vertices[i];
            Vector2D edgeB = B.vertices[j+1]-B.vertices[j];
            cross=edgeA.Cross(edgeB);
        }
        if(j==m||(i<n&&Sign(cross)>=0)){
            ++i;
        }
        if(i==n||(j<m&&Sign(cross)<=0)){
            ++j;
        }
    }

    //移除多补的一个末尾顶点
    A.vertices.pop_back();
    B.vertices.pop_back();

    return result;
}


int n1 = 0, n2 = 0, m = 0;
Polygon polygon1;
Polygon polygon2;
vector<Vector2D> testCases;

struct Projection {
    double min, max;
};

Projection ProjectPolygon(const Polygon& poly, const Vector2D& axis)
{
    double minProj = poly.vertices[0].Dot(axis);
    double maxProj = minProj;

    for (size_t i = 1; i < poly.vertices.size(); ++i) {
        double proj = poly.vertices[i].Dot(axis);
        if (proj < minProj) {
            minProj = proj;
        }
        if (proj > maxProj) {
            maxProj = proj;
        }
    }
    return {minProj, maxProj};
}
// Vector2D GenSolution(const Vector2D& vec)
// {
//     Polygon polyB = polygon2;
//     polyB.MoveByVec(vec);

//     double minOverlap = std::numeric_limits<double>::infinity();
//     Vector2D smallestAxis;

//     const Polygon* polygons[2] = {&polygon1, &polyB};

//     for (int i = 0; i < 2; ++i) {
//         const Polygon& currentPoly = *polygons[i];
//         for (size_t j = 0; j < currentPoly.vertices.size(); ++j) {
//             Vector2D p1 = currentPoly.vertices[j];
//             Vector2D p2 = currentPoly.vertices[(j + 1) % currentPoly.vertices.size()];
//             Vector2D edge = p2 - p1;

//             Vector2D axis = edge.Perp().Normalize();

//             Projection projA = ProjectPolygon(polygon1, axis);
//             Projection projB = ProjectPolygon(polyB, axis);

//             if (projA.max <= projB.min || projB.max <= projA.min) {
//                 return {0.0, 0.0};
//             }

//             double overlap = std::min(projA.max, projB.max) - std::max(projA.min, projB.min);

//             if (overlap < minOverlap) {
//                 minOverlap = overlap;
//                 smallestAxis = axis;
//             }
//         }
//     }
//     Vector2D centerA = polygon1.GetCenter();
//     Vector2D centerB = polyB.GetCenter();
//     Vector2D dir = centerB - centerA;

//     if (smallestAxis.Dot(dir) < 0.0) {
//         smallestAxis = smallestAxis * -1.0;
//     }

//     return {smallestAxis * minOverlap};
// }
Vector2D GenSolution(const Vector2D& vec){

}

// 选手在规定的时间内进行预处理，完成后返回OK
void PreProcess()
{
    // player can perform preprocessing here
}

int main()
{
    // =============== 1. read polygons ===================
    cin >> n1 >> n2;

    if (n1 <= 0 || n2 <= 0) {
        cerr << "Input data error: the number of vertices of both polygons should be greater than 2" << endl;
        return 1;
    }

    polygon1.vertices.resize(n1);
    for (int i = 0; i < n1; ++i) {
        cin >> polygon1.vertices[i].x >> polygon1.vertices[i].y;
    }

    polygon2.vertices.resize(n2);
    for (int i = 0; i < n2; ++i) {
        cin >> polygon2.vertices[i].x >> polygon2.vertices[i].y;
    }

    // wait for OK to ensure all polygon data is received
    string okResp;
    cin >> okResp;
    if (okResp != "OK") {
        cerr << "Input data error: waiting for OK after obtaining polygons but I get " << okResp << endl;
        return 0;
    }

    // ============== 2. preprocess ===================
    PreProcess();
    // send OK after finishing preprocessing
    cout << "OK" << endl;
    cout.flush();

    // ============== 3. read test data ===================
    cin >> m;
    testCases.resize(m);

    for (int i = 0; i < m; ++i) {
        cin >> testCases[i].x >> testCases[i].y;
    }

    // wait for OK to ensure all test cases are received
    cin >> okResp;
    if (okResp != "OK") {
        cerr << "Input data error: waiting for OK after that I have received all test points but I get " << okResp
             << endl;
        return 0;
    }

    // ================ 4. solve and output results ===================
    for (int i = 0; i < m; ++i) {
        const Vector2D& res = GenSolution(testCases[i]);
        cout << fixed << std::setprecision(5) << res.x << " " << res.y << endl;
        cout.flush();
    }

    // send OK after outputting all answer
    cout << "OK" << endl;
    cout.flush();

    return 0;
}