/*H***********************************************************
 *
    Header file for vector
 *
 *H*/

#ifndef VECTOR_H_
#define VECTOR_H_

typedef struct {double x, y;} VecR2;
typedef struct {double x, y, z;} VecR3;

#define Sqr(x) ((x) * (x))
#define Cub(x) ((x) * (x) * (x))

// 2D
#define VSet2(v, sx, sy)             \
(v).x = sx,                          \
(v).y = sy
#define VAdd2(v1, v2, v3)            \
(v1).x = (v2).x + (v3).x,            \
(v1).y = (v2).y + (v3).y
#define VSub2(v1, v2, v3)            \
(v1).x = (v2).x - (v3).x,            \
(v1).y = (v2).y - (v3).y
#define VMul2(v1, v2, v3)            \
(v1).x = (v2).x * (v3).x,            \
(v1).y = (v2).y * (v3).y
#define VDiv2(v1, v2, v3)            \
(v1).x = (v2).x / (v3).x,            \
(v1).y = (v2).y / (v3).y
#define VDot2(v1, v2)                \
((v1).x * (v2).x + (v1).y * (v2).y)
#define VSAdd2(v1, v2, s, v3)        \
(v1).x = (v2).x + (s)* (v3).x,	     \
(v1).y = (v2).y + (s)* (v3).y
#define VSCopy2(v1, s1, v2)          \
(v1).x = (s1)* (v2).x,				 \
(v1).y = (s1)* (v2).y

// 3D
#define VSet3(v, sx, sy, sz)               \
v.x = sx,                                  \
v.y = sy,							       \
v.z = sz
#define VAdd3(v1, v2, v3)                  \
(v1).x = (v2).x + (v3).x,                  \
(v1).y = (v2).y + (v3).y,		           \
(v1).z = (v2).z + (v3).z
#define VSub3(v1, v2, v3)                  \
(v1).x = (v2).x - (v3).x,                  \
(v1).y = (v2).y - (v3).y,		           \
(v1).z = (v2).z - (v3).z
#define VMul3(v1, v2, v3)                  \
(v1).x = (v2).x * (v3).x,                  \
(v1).y = (v2).y * (v3).y,				   \
(v1).z = (v2).z * (v3).z
#define VDiv3(v1, v2, v3)                  \
(v1).x = (v2).x / (v3).x,                  \
(v1).y = (v2).y / (v3).y,				   \
(v1).z = (v2).z / (v3).z,
#define VDot3(v1, v2)                      \
(v1.x * v2.x + v1.y * v2.y + v1.z * v2.z)
#define VSAdd3(v1, v2, s, v3)              \
(v1).x = (v2).x + (s)* (v3).x,			   \
(v1).y = (v2).y + (s)* (v3).y,			   \
(v1).z = (v2).z + (s)* (v3).z
#define VSCopy3(v1, s1, v2)                \
(v1).x = (s1)* (v2).x,					   \
(v1).y = (s1)* (v2).y,					   \
(v1).z = (s1)* (v2).z

#endif
