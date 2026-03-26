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

// 【精度极限】将容差收紧至 1e-9，并全盘使用 long double
const long double EPS = 1e-9;

// 精度控制函数
inline int Sign(long double x){
    if(x > EPS) return 1;
    if(x < -EPS) return -1;
    return 0;
}

// 实现二维向量的基本操作
struct Vector2D {
    long double x, y;

    Vector2D(long double x = 0, long double y = 0) : x(x), y(y) {}
    
    // 【性能与精度】传值优化，配合 80位 long double 寄存器
    Vector2D operator-(Vector2D other) const { return Vector2D(x - other.x, y - other.y); }
    Vector2D operator+(Vector2D other) const { return Vector2D(x + other.x, y + other.y); }
    Vector2D operator*(long double scalar) const { return Vector2D(x * scalar, y * scalar); }
    
    long double Dot(Vector2D other) const { return x * other.x + y * other.y; }
    long double Cross(Vector2D other) const { return x * other.y - y * other.x; }
    long double Length() const { return std::sqrt(x * x + y * y); }
    
    Vector2D Normalize() const {
        long double len = Length();
        if (len == 0) return *this;
        return Vector2D(x / len, y / len);
    }
    Vector2D Perp() const { return Vector2D(-y, x); }
};

inline long double CrossProduct(Vector2D a, Vector2D b, Vector2D c) {
    return (b - a).Cross(c - a);
}

// 点到线段最短距离
struct DistanceInfo{
    long double distanceSq; 
    Vector2D pushVector;
};

// 【精度极限】引入端点精确钳位 (Exact Clamping)
DistanceInfo PointToSegmentDistance(Vector2D p, Vector2D a, Vector2D b){
    Vector2D ab = b - a;
    Vector2D ap = p - a;
    long double ab_len_sq = ab.Dot(ab);
    
    if (Sign(ab_len_sq) == 0) {
        Vector2D p_to_a = a - p;
        return {p_to_a.Dot(p_to_a), p_to_a};
    }
    
    long double t = ap.Dot(ab) / ab_len_sq;
    
    // 如果投影落在端点 a 外侧，直接返回 a-p，避免浮点乘加误差
    if (t <= 0.0) {
        Vector2D p_to_a = a - p;
        return {p_to_a.Dot(p_to_a), p_to_a};
    }
    // 如果投影落在端点 b 外侧，直接返回 b-p
    if (t >= 1.0) {
        Vector2D p_to_b = b - p;
        return {p_to_b.Dot(p_to_b), p_to_b};
    }
    
    // 严格投影在线段内部
    Vector2D closestPoint = a + ab * t;
    Vector2D pushVector = closestPoint - p;
    return {pushVector.Dot(pushVector), pushVector};
}

struct Polygon {
    std::vector<Vector2D> vertices;
    long double minX = 0, minY = 0, maxX = 0, maxY = 0; 

    Polygon() = default;
    Polygon(std::vector<Vector2D> vts) : vertices(vts) {
        UpdateAABB();
    }

    void UpdateAABB() {
        if (vertices.empty()) return;
        minX = maxX = vertices[0].x;
        minY = maxY = vertices[0].y;
        for (const auto& v : vertices) {
            minX = std::min(minX, v.x);  maxX = std::max(maxX, v.x);
            minY = std::min(minY, v.y);  maxY = std::max(maxY, v.y);
        }
    }

    // 【鲁棒性优化】凸多边形包含测试，不依赖射线法，不受精度抖动干扰
    bool ContainsConvex(Vector2D p) const {
        if (p.x < minX - EPS || p.x > maxX + EPS || 
            p.y < minY - EPS || p.y > maxY + EPS) {
            return false;
        }

        int n = vertices.size();
        if (n < 3) return false;
        for (int i = 0; i < n; ++i) {
            // 只要点在任何一条边的右侧，即在多边形外
            if (Sign(CrossProduct(vertices[i], vertices[(i + 1) % n], p)) < 0) {
                return false; 
            }
        }
        return true; 
    }
};

struct Segment {
    Vector2D a, b;
    Segment(Vector2D a, Vector2D b) : a(a), b(b) {}
};

// 获得对称多边形
Polygon GetNativePolygon(const Polygon& poly){
    Polygon negPoly;
    negPoly.vertices.reserve(poly.vertices.size());
    for(const auto& v: poly.vertices){
        negPoly.vertices.emplace_back(-v.x, -v.y);
    }
    negPoly.UpdateAABB(); 
    return negPoly;
}

bool IsPointInTiangle(Vector2D p, Vector2D a, Vector2D b, Vector2D c){
    long double cp1 = CrossProduct(a, b, p);
    long double cp2 = CrossProduct(b, c, p);
    long double cp3 = CrossProduct(c, a, p);
    return Sign(cp1) >= 0 && Sign(cp2) >= 0 && Sign(cp3) >= 0;
}

// 耳切法实现多边形分解
std::vector<Polygon> PartPolygonToTriangle(const Polygon& poly) {
    std::vector<Polygon> triangles;
    std::vector<Vector2D> V = poly.vertices;
    int n = V.size();

    if (n <= 3) {
        triangles.emplace_back(poly);
        return triangles;
    }

    std::vector<bool> isReflex(n, false);
    for (int i = 0; i < n; i++) {
        int prev = (i - 1 + n) % n;
        int next = (i + 1) % n;
        if (Sign(CrossProduct(V[prev], V[i], V[next])) <= 0) {
            isReflex[i] = true;
        }
    }

    std::vector<int> nextIdx(n), prevIdx(n);
    for (int i = 0; i < n; i++) {
        nextIdx[i] = (i + 1) % n;
        prevIdx[i] = (i - 1 + n) % n;
    }

    int remain = n;
    int curr = 0;
    int maxAttempts = 2 * n; 

    while (remain > 3 && maxAttempts > 0) {
        maxAttempts--;
        int prev = prevIdx[curr];
        int next = nextIdx[curr];

        if (Sign(CrossProduct(V[prev], V[curr], V[next])) > 0) {
            bool isEar = true;
            int testIdx = nextIdx[next];
            while (testIdx != prev) {
                if (isReflex[testIdx]) {
                    if (IsPointInTiangle(V[testIdx], V[prev], V[curr], V[next])) {
                        isEar = false;
                        break;
                    }
                }
                testIdx = nextIdx[testIdx];
            }

            if (isEar) {
                triangles.emplace_back(Polygon(std::vector<Vector2D>{V[prev], V[curr], V[next]}));
                nextIdx[prev] = next;
                prevIdx[next] = prev;
                remain--;
                int prevPrev = prevIdx[prev];
                isReflex[prev] = Sign(CrossProduct(V[prevPrev], V[prev], V[next])) <= 0;
                int nextNext = nextIdx[next];
                isReflex[next] = Sign(CrossProduct(V[prev], V[next], V[nextNext])) <= 0;
                curr = prev;
                maxAttempts = 2 * remain; 
                continue;
            }
        }
        curr = nextIdx[curr];
    }

    if (remain == 3) {
        int v1 = curr;
        int v2 = nextIdx[v1];
        int v3 = nextIdx[v2];
        triangles.emplace_back(Polygon(std::vector<Vector2D>{V[v1], V[v2], V[v3]}));
    }

    return triangles;
}

void Reorder(std::vector<Vector2D>& points){
    int minidx = 0;
    int n = points.size();
    for(int i = 1; i < n; i++){
        if(points[i].y < points[minidx].y || (Sign(points[i].y - points[minidx].y) == 0 && points[i].x < points[minidx].x)){
            minidx = i;
        }
    }
    std::rotate(points.begin(), points.begin() + minidx, points.end());
}

Polygon MinkowskiSum(Polygon A, Polygon B){
    Reorder(A.vertices);
    Reorder(B.vertices);
    
    A.vertices.push_back(A.vertices[0]);
    B.vertices.push_back(B.vertices[0]);

    Polygon result;
    int i = 0, j = 0;
    int n = A.vertices.size() - 1;
    int m = B.vertices.size() - 1;
    
    result.vertices.reserve(n + m);

    while(i < n || j < m){
        result.vertices.push_back(A.vertices[i] + B.vertices[j]);
        long double cross = 0;
        if(i < n && j < m){
            Vector2D edgeA = A.vertices[i+1] - A.vertices[i];
            Vector2D edgeB = B.vertices[j+1] - B.vertices[j];
            cross = edgeA.Cross(edgeB);
        }
        if(j == m || (i < n && Sign(cross) >= 0)){ ++i; }
        if(i == n || (j < m && Sign(cross) <= 0)){ ++j; }
    }

    A.vertices.pop_back();
    B.vertices.pop_back();
    result.UpdateAABB(); 
    return result;
}

// ====================== 全局环境 ======================
int n1 = 0, n2 = 0, m = 0;
Polygon polygon1;
Polygon polygon2;
vector<Vector2D> testCases;

std::vector<Polygon> nfpFragments;
std::vector<Segment> unionBoundarySegments;

long double globalNfpMinX = std::numeric_limits<long double>::max();
long double globalNfpMaxX = std::numeric_limits<long double>::lowest();
long double globalNfpMinY = std::numeric_limits<long double>::max();
long double globalNfpMaxY = std::numeric_limits<long double>::lowest();

bool IsInsideUnionNFP(Vector2D p) {
    if (p.x < globalNfpMinX - EPS || p.x > globalNfpMaxX + EPS || 
        p.y < globalNfpMinY - EPS || p.y > globalNfpMaxY + EPS) {
        return false;
    }
    for (const auto& poly : nfpFragments) {
        if (poly.ContainsConvex(p)) return true;
    }
    return false;
}

// ====================== 核心查询流程 ======================
Vector2D GenSolution(Vector2D vec) {
    if (unionBoundarySegments.empty()) return Vector2D(0.0, 0.0);

    // 1. 重叠检查
    if (!IsInsideUnionNFP(vec)) return Vector2D(0.0, 0.0);

    // 2. 试探验证与最短距离查找
    long double minDistanceSq = std::numeric_limits<long double>::infinity();
    Vector2D bestPushVector(0, 0);
    
    // 鲁棒兜底
    long double fallbackDistSq = std::numeric_limits<long double>::infinity();
    Vector2D fallbackPush(0, 0);

    for (const auto& seg : unionBoundarySegments) {
        DistanceInfo info = PointToSegmentDistance(vec, seg.a, seg.b);
        long double distSq = info.distanceSq;
        Vector2D push_vec = info.pushVector;

        if (distSq < fallbackDistSq) {
            fallbackDistSq = distSq;
            fallbackPush = push_vec;
        }

        if (distSq >= minDistanceSq) continue;

        if (distSq < EPS * EPS) {
            minDistanceSq = 0;
            bestPushVector = Vector2D(0.0, 0.0);
            continue;
        }

        // 【极致精度】缩小步长至 1e-7 进行外法线试探
        Vector2D normal = push_vec.Normalize();
        Vector2D testPoint = vec + push_vec + normal * 1e-7; 
        
        bool isTargetValid = true;
        for (const auto& nfp : nfpFragments) {
            if (nfp.ContainsConvex(testPoint)) {
                isTargetValid = false;
                break;
            }
        }
        
        if (isTargetValid) {
            minDistanceSq = distSq;
            bestPushVector = push_vec;
        }
    }

    if (minDistanceSq == std::numeric_limits<long double>::infinity()) {
        return fallbackPush;
    }

    return bestPushVector;
}

void PreProcess()
{
    Polygon negB = GetNativePolygon(polygon2);
    std::vector<Polygon> partsA = PartPolygonToTriangle(polygon1);
    std::vector<Polygon> partsNegB = PartPolygonToTriangle(negB);

    for(const auto& partA : partsA){
        for(const auto& partB : partsNegB){
            nfpFragments.push_back(MinkowskiSum(partA, partB));
        }
    }

    for (const auto& poly : nfpFragments) {
        int n = poly.vertices.size();
        for (int i = 0; i < n; ++i) {
            Vector2D p1 = poly.vertices[i];
            Vector2D p2 = poly.vertices[(i + 1) % n];
            unionBoundarySegments.emplace_back(p1, p2);

            globalNfpMinX = std::min({globalNfpMinX, p1.x, p2.x});
            globalNfpMaxX = std::max({globalNfpMaxX, p1.x, p2.x});
            globalNfpMinY = std::min({globalNfpMinY, p1.y, p2.y});
            globalNfpMaxY = std::max({globalNfpMaxY, p1.y, p2.y});
        }
    }
}

int main()
{
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    if (!(cin >> n1 >> n2)) return 0;

    polygon1.vertices.resize(n1);
    for (int i = 0; i < n1; ++i) {
        double tx, ty; cin >> tx >> ty;
        polygon1.vertices[i] = Vector2D(tx, ty);
    }
    polygon1.UpdateAABB(); 

    polygon2.vertices.resize(n2);
    for (int i = 0; i < n2; ++i) {
        double tx, ty; cin >> tx >> ty;
        polygon2.vertices[i] = Vector2D(tx, ty);
    }
    polygon2.UpdateAABB(); 

    string okResp;
    cin >> okResp;
    if (okResp != "OK") return 0;

    PreProcess();
    cout << "OK\n"; cout.flush();

    cin >> m;
    testCases.resize(m);
    for (int i = 0; i < m; ++i) {
        double tx, ty; cin >> tx >> ty;
        testCases[i] = Vector2D(tx, ty);
    }

    cin >> okResp;
    if (okResp != "OK") return 0;

    cout << m << "\n";
    for (int i = 0; i < m; ++i) {
        Vector2D res = GenSolution(testCases[i]);
        // 输出 6 位小数以匹配精度需求
        cout << fixed << std::setprecision(6) << (double)res.x << " " << (double)res.y << "\n";
    }

    cout << "OK\n"; cout.flush();
    return 0;
}