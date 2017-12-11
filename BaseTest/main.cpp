#include <iostream>
#include <vector>
#include <map>
#include <unordered_map>
#include <cmath>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

#include <glm.hpp>

typedef OpenMesh::TriMesh_ArrayKernelT<> MyMesh;

typedef std::vector<std::pair<MyMesh::EdgeHandle, bool>> DeletingEdgeArray;
typedef std::vector<std::pair<MyMesh::EdgeHandle, bool>> UnDeterminedEdgeArray;
typedef std::vector<std::pair<MyMesh::EdgeHandle, bool>> DeletedEdgeArray;

#define FLOAT_MAX 1000000000000.0f

#define DELETING 0
#define DELETED 1
#define UNDETERMINED 2
#define BORDER 3
#define BORDER_CALCED 4

//struct HashFunc
//{
//	std::size_t operator()(const MyMesh::EdgeHandle &key) const
//	{
//		return key.idx();
//	}
//};
//
//struct EqualKey
//{
//	bool operator() (const MyMesh::EdgeHandle &lhs, const MyMesh::EdgeHandle & rhs)
//	{
//		return lhs.idx() == rhs.idx();
//	}
//};
//
//std::unordered_map<MyMesh::EdgeHandle, int, HashFunc, EqualKey> mp;

typedef struct Point_
{
	Point_(float _x = 0.0f, float _y = 0.0f, float _z = 0.0f) :
		x(_x),
		y(_y),
		z(_z)
	{
	};

	Point_(const Point_ & p)
	{
		x = p.x;
		y = p.y;
		z = p.z;
	}

	bool operator< (const Point_ & rhs)
	{
		return this->x < rhs.x
			|| (this->x == rhs.x) && (this->y < rhs.y)
			|| (this->x == rhs.x && this->y == rhs.y) && (this->z < rhs.z);
	}

	bool operator==(const Point_ & rhs)
	{
		return this->x == rhs.x && this->y == rhs.y && this->z == rhs.z;
	}

	bool operator() (const Point_ & lhs, const Point_ & rhs)
	{
		return lhs.x < rhs.x
			|| (lhs.x == rhs.x) && (lhs.y < rhs.y)
			|| (lhs.x == rhs.x && lhs.y == rhs.y) && (lhs.z < rhs.z);
	}

	float x;
	float y;
	float z;
} Point;

typedef Point BorderPoint;

typedef struct Region_
{
	Region_() : isDeleted(false) {};

	std::vector<MyMesh::FaceHandle> faces;
	typedef int blsIdx;
	std::set<blsIdx> blss;
	bool isDeleted;
} Region;

typedef struct BorderLineSegment_
{
	BorderLineSegment_(Point& _A, const Point& _B)
	{
		if (!(_A < _B))
		{
			A = _B;
			B = _A;
		}
		else
		{
			A = _A;
			B = _B;
		}

		length = glm::length(glm::vec3(A.x - B.x, A.y - B.y, A.z - B.z));
		needsDeleted = false;
	}

	Point A;
	Point B;

	float length;
	std::set<int> RegionIDs;

	bool needsDeleted;
} BorderLineSegment;

typedef int BLS_id; // Border Line Segments idx
typedef struct GraphNode_
{
	GraphNode_() : isVisited(false) {};

	std::vector<std::pair<BorderPoint, BLS_id>> toPoints;
	bool isVisited;
} GraphNode;

//生成边界点构成的最大连通区域
typedef struct ConnectRegion_
{
	ConnectRegion_() : isDeleted(false), meanLength(FLOAT_MAX) {};

	std::vector<BorderPoint> points;
	typedef int blsIdx;
	std::set<blsIdx> blss;
	float meanLength;
	bool isDeleted;

} ConnectRegion;


void RegionSpread(MyMesh & mesh, OpenMesh::EPropHandleT<int>& edgeStatus, MyMesh::FaceHandle &fh, std::map<MyMesh::FaceHandle, bool> &mp);
void GraphBFS(std::map < BorderPoint, GraphNode, BorderPoint>& mp, BorderPoint & startPoint, ConnectRegion & newRegion);
void DeletedRegionAnalysis(std::vector<ConnectRegion> & cr, std::map < BorderPoint, GraphNode, BorderPoint>& mp, std::vector<BorderLineSegment> &borderLineSegments);

void DeleteMyMesh(MyMesh & mesh,
	std::vector<ConnectRegion> & cr,
	std::vector<BorderLineSegment> &borderLineSegments,
	std::vector<Region>& regionPs);

int main(int argc, char **argv)
{
	MyMesh mesh;

	//control parameter
	if (argc != 2)
	{
		std::cerr << "Please input the mesh file name! " << std::endl;
		system("pause");
		return 1;
	}

	//read mesh
	if (!OpenMesh::IO::read_mesh(mesh, argv[1]))
	{
		std::cerr << "Error: Cannot read mesh from " << argv[1] << std::endl;
		system("pause");
		return 1;
	}

	OpenMesh::IO::Options opt;
	//add property
	//normal property
	if (!opt.check(OpenMesh::IO::Options::FaceNormal))
	{
		mesh.request_face_normals();
		mesh.update_normals();
	}

	//dihedral angle
	OpenMesh::EPropHandleT<float> dihedralAngle;
	mesh.add_property(dihedralAngle);

	//Edge status 
	OpenMesh::EPropHandleT<int> status;
	mesh.add_property(status);

	OpenMesh::EPropHandleT<std::set<int>> regionIDs;
	mesh.add_property(regionIDs);

	//平面表面区域边界提取及分割
	//计算每一条边相邻面的二面角
	DeletingEdgeArray deletingEdges;
	UnDeterminedEdgeArray undeterminedEdges;

	float fi = 0.1f;
	for (MyMesh::EdgeIter e_it = mesh.edges_begin();
		e_it != mesh.edges_end(); ++e_it)
	{
		float f = mesh.calc_dihedral_angle_fast(*e_it);
		float degree = glm::degrees(f);
		mesh.property(dihedralAngle, *e_it) = degree;
		if (std::abs(degree) < fi)
		{
			mesh.property(status, *e_it) = DELETING;
			deletingEdges.push_back(std::make_pair(*e_it, true));
		}
		else
		{
			mesh.property(status, *e_it) = UNDETERMINED;
			undeterminedEdges.push_back(std::make_pair(*e_it, true));
		}
	}

	//生成平面区域
	std::vector<Region> regionPs;
	for (int i = 0; i < deletingEdges.size(); ++i)
	{
		std::map<MyMesh::FaceHandle, bool> mp;
		if (mesh.property(status, deletingEdges[i].first) == DELETING) //边状态为待删除
		{
			mp.clear();
			Region newRegion;
			auto adjacentfh = mesh.face_handle(mesh.halfedge_handle(deletingEdges[i].first, 0));
			RegionSpread(mesh, status, adjacentfh, mp);
			for (const auto & i : mp)
			{
				newRegion.faces.push_back(i.first);
			}
			regionPs.push_back(newRegion);
		}
	}

	int borderSum = 0;
	for (int i = 0; i < regionPs.size(); ++i)
	{
		for (int j = 0; j < regionPs[i].faces.size(); ++j)
		{
			MyMesh::FaceEdgeIter fe_iter = mesh.fe_iter(regionPs[i].faces[j]);
			while (fe_iter.is_valid())
			{
				if (mesh.property(status, *fe_iter) == UNDETERMINED)
				{
					mesh.property(status, *fe_iter) = BORDER;
					mesh.property(regionIDs, *fe_iter).insert(i);
					++borderSum;
				}
				else if (mesh.property(status, *fe_iter) == BORDER)
				{
					mesh.property(regionIDs, *fe_iter).insert(i);
				}

				++fe_iter;
			}
		}
	}

	std::cout << "Region Sum: " << regionPs.size() << std::endl;
	std::cout << "Border Sum: " << borderSum << std::endl;

	std::vector<BorderLineSegment> borderLineSegments;
	for (int i = 0; i < regionPs.size(); ++i)
	{
		for (int j = 0; j < regionPs[i].faces.size(); ++j)
		{
			MyMesh::FaceEdgeIter fe_iter = mesh.fe_iter(regionPs[i].faces[j]);
			while (fe_iter.is_valid())
			{
				if (mesh.property(status, *fe_iter) == BORDER)
				{
					MyMesh::Point pA = 
						mesh.point(mesh.from_vertex_handle(mesh.halfedge_handle(*fe_iter, 0)));
					MyMesh::Point pB = 
						mesh.point(mesh.to_vertex_handle(mesh.halfedge_handle(*fe_iter, 0)));
					std::cout << "From : " << pA[0] << " " << pA[1] << " " << pA[2]
						<< " To : " << pB[0] << " " << pB[1] << " " << pB[2] << std::endl;

					BorderLineSegment bls(BorderLineSegment(
						Point(pA[0], pA[1], pA[2]),
						Point(pB[0], pB[1], pB[2])));

					for (const auto & i : mesh.property(regionIDs, *fe_iter))
					{
						bls.RegionIDs.insert(i);
					}

					borderLineSegments.push_back(bls);
					mesh.property(status, *fe_iter) = BORDER_CALCED;
				}

				++fe_iter;
			}
		}
	}
	
	std::map < BorderPoint, GraphNode, BorderPoint> borderPointAdjGraph;
	for (int i = 0; i < borderLineSegments.size(); ++i)
	{
		borderPointAdjGraph[borderLineSegments[i].A].toPoints.push_back({ borderLineSegments[i].B, i });
		borderPointAdjGraph[borderLineSegments[i].B].toPoints.push_back({ borderLineSegments[i].A, i });

		for (const auto & j : borderLineSegments[i].RegionIDs)
		{
			regionPs[j].blss.insert(i);
		}
	}

	std::vector<ConnectRegion> connectRegions;
	for (auto & bp : borderPointAdjGraph)
	{
		if (bp.second.isVisited) continue;
		ConnectRegion newRegion;
		BorderPoint startPoint(bp.first);
		GraphBFS(borderPointAdjGraph, startPoint, newRegion);
		connectRegions.push_back(newRegion);
	}
	
	DeletedRegionAnalysis(connectRegions, borderPointAdjGraph, borderLineSegments);

	DeleteMyMesh(mesh, connectRegions, borderLineSegments, regionPs);

	if (!OpenMesh::IO::write_mesh(mesh, "result.off"))
	{
		std::cerr << "write error \n";
		return 1;
	}

	//delete faces edges vertices
	/*mesh.request_face_status();
	mesh.request_edge_status();
	mesh.request_vertex_status();

	mesh.delete_face(fh, true);
	mesh.delete_vertex(vh, true);
	mesh.delete_edge(eh, true);

	mesh.garbage_collection();*/



	/*std::cout << "Before: " << std::endl;
	for (MyMesh::EdgeIter e_it = mesh.edges_begin();
	e_it != mesh.edges_end(); ++e_it)
	{
	std::cout << "Edge #" << *e_it << " : " << mesh.property(dihedralAngle, *e_it) << std::endl;
	}

	std::cout << "After : " << std::endl;
	for (MyMesh::EdgeIter e_it = mesh.edges_begin();
	e_it != mesh.edges_end(); ++e_it)
	{
	std::cout << "Edge #" << *e_it << " : " << mesh.property(dihedralAngle, *e_it) << std::endl;

	auto fh1 = mesh.face_handle(mesh.halfedge_handle(*e_it, 0));
	auto fh2 = mesh.opposite_face_handle(mesh.halfedge_handle(*e_it, 0));
	}*/


	system("pause");
	return 0;
}

void RegionSpread(MyMesh & mesh, OpenMesh::EPropHandleT<int>& edgeStatus, MyMesh::FaceHandle &fh, std::map<MyMesh::FaceHandle, bool> &mp)
{
	if (mp.count(fh) > 0) return;

	mp.insert(std::make_pair(fh, true));

	MyMesh::FaceHalfedgeIter fhe_it = mesh.fh_iter(fh);
	
	while (fhe_it.is_valid())
	{
		MyMesh::HalfedgeHandle heh = *fhe_it;
		MyMesh::EdgeHandle eh = mesh.edge_handle(heh);
		if (mesh.property(edgeStatus, eh) == DELETING)
		{
			mesh.property(edgeStatus, eh) = DELETED;
			MyMesh::FaceHandle fhNew = mesh.opposite_face_handle(heh);
			RegionSpread(mesh, edgeStatus, fhNew, mp);
		}
		++fhe_it;
	}
}

void GraphBFS(std::map < BorderPoint, GraphNode, BorderPoint>& mp, BorderPoint & startPoint, ConnectRegion & newRegion)
{
	if (mp[startPoint].isVisited) return;

	mp[startPoint].isVisited = true;
	newRegion.points.push_back(startPoint);

	for (const auto & i : mp[startPoint].toPoints)
	{
		newRegion.blss.insert(i.second);
	}

	for (auto & adjPoint : mp[startPoint].toPoints)
	{
		GraphBFS(mp, adjPoint.first, newRegion);
	}
}

void DeletedRegionAnalysis(std::vector<ConnectRegion> & cr, std::map < BorderPoint, GraphNode, BorderPoint>& mp, std::vector<BorderLineSegment> &borderLineSegments)
{
	for (auto & r : cr)
	{
		float meanLength = 0.0f;
		int n = 0;
		for (auto iter = r.blss.begin(); iter != r.blss.end(); ++iter, ++n)
		{
			meanLength += borderLineSegments[*iter].length;
		}
		meanLength /= n;
		r.meanLength = meanLength;
	}

	sort(cr.begin(), cr.end(), [](const ConnectRegion& lhs, const ConnectRegion& rhs){
		return lhs.meanLength < rhs.meanLength;
	});

	cr[0].isDeleted = true;
	cr[1].isDeleted = true;
}

void DeleteMyMesh(MyMesh & mesh, 
	std::vector<ConnectRegion> & cr, 
	std::vector<BorderLineSegment> &borderLineSegments, 
	std::vector<Region>& regionPs)
{
	for (const auto & crAuto : cr)
	{
		if (crAuto.isDeleted)
		{
			for (const auto & blssAuto : crAuto.blss)
			{
				borderLineSegments[blssAuto].needsDeleted = true;
			}
		}
	}

	for (auto & rPsAuto : regionPs)
	{
		bool isDeleted = false;
		auto iter = rPsAuto.blss.begin();
		for (; iter != rPsAuto.blss.end(); ++iter)
		{
			if (!borderLineSegments[*iter].needsDeleted) break;
		}
		if (iter == rPsAuto.blss.end()) isDeleted = true;
		rPsAuto.isDeleted = isDeleted;
	}

	mesh.request_face_status();
	mesh.request_edge_status();
	mesh.request_vertex_status();

	for (int i = 0; i < regionPs.size(); ++i)
	{
		if (regionPs[i].isDeleted)
		{
			for (auto & fh : regionPs[i].faces)
				mesh.delete_face(fh, true);
		}
	}
	mesh.garbage_collection();
}