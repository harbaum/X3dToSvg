#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, os, math, copy
import xml.etree.ElementTree as et

import pysvg.structure
import pysvg.builders
import pysvg.text

# ============================== depth buffer ================================

WIDTH=1000
HEIGHT=1000
zbuf = [[None for i in range(WIDTH)] for j in range(HEIGHT)]
coverbuf = { }

def drawHLine(x0,x1,y,z, f):
    if x0 > x1:
        t = x0
        x0 = x1
        x1 = t

    for x in range(x0,x1+1):
        if zbuf[y][x] != None:
            if z > zbuf[y][x][0]:
                if not f["id"] in coverbuf:
                    coverbuf[f["id"]] = []
                if not zbuf[y][x][1]["id"] in coverbuf[f["id"]]:
                    coverbuf[f["id"]].append(zbuf[y][x][1]["id"])
            else:
                if not zbuf[y][x][1]["id"] in coverbuf:
                    coverbuf[zbuf[y][x][1]["id"]] = []
                if not f["id"] in coverbuf[zbuf[y][x][1]["id"]]:
                    coverbuf[zbuf[y][x][1]["id"]].append(f["id"])
        
        if zbuf[y][x] == None or z > zbuf[y][x][0]:
            zbuf[y][x] = [ z, f ]

# http://www.sunshine2k.de/coding/java/TriangleRasterization/TriangleRasterization.html
def bottomTriangle(v1, v2, v3, f):
    # v2[1] == v3[1]
    # print(v1,v2,v3)

    # just one line?
    if v2[1] == v1[1]:
        drawHLine(v2[0], v3[0], v1[1], v2[2], f)
        return
    
    invslope1 = (v2[0] - v1[0]) / (v2[1] - v1[1]);
    invslope2 = (v3[0] - v1[0]) / (v3[1] - v1[1]);
    zstep     = (v2[2] - v1[2]) / (v2[1] - v1[1]);
    
    curx1 = v1[0];
    curx2 = v1[0];
    curz  = v1[2];

    scanlineY = int(v1[1])
    while scanlineY <= v2[1]:
        drawHLine(int(curx1), int(curx2), scanlineY, curz ,f);
        curx1 = curx1 + invslope1;
        curx2 = curx2 + invslope2;
        curz = curz + zstep;
        scanlineY = scanlineY + 1

def topTriangle(v1, v2, v3, f):
    # v1[1] == v2[1]
    
    # just one line?
    if v2[1] == v3[1]:
        drawHLine(v1[0], v2[0], v1[1], v2[2], f)
        return
    
    invslope1 = (v3[0] - v1[0]) / (v3[1] - v1[1]);
    invslope2 = (v3[0] - v2[0]) / (v3[1] - v2[1]);
    zstep     = (v3[2] - v1[2]) / (v3[1] - v1[1]);

    curx1 = v3[0];
    curx2 = v3[0];
    curz  = v3[2];

    scanlineY = int(v3[1])
    while scanlineY > v1[1]:
        drawHLine(int(curx1), int(curx2), scanlineY, curz, f);
        curx1 = curx1 - invslope1;
        curx2 = curx2 - invslope2;
        curz = curz - zstep;        
        scanlineY = scanlineY - 1

def ysort(c):
    return c[1]
        
def drawTriangle(t, f):
    # at first sort the three vertices by y-coordinate ascending so v[1] is the topmost vertice 
    # sortVerticesAscendingByY();
    v = sorted(t, key=ysort)
    
    # here we know that v1[1] <= v2[1] <= v3[1]
    # check for trivial case of bottom-flat triangle
    if v[1][1] == v[2][1]:
        bottomTriangle(v[0], v[1], v[2], f);

    # check for trivial case of top-flat triangle
    elif v[0][1] == v[1][1]:
        topTriangle(v[0], v[1], v[2], f);

    # general case - split the triangle in a topflat and bottom-flat one
    else:
        v3 = [ int((v[0][0] + ((v[1][1] - v[0][1]) / (v[2][1] - v[0][1])) * (v[2][0] - v[0][0]))), int(v[1][1]),
               int((v[0][2] + ((v[1][1] - v[0][1]) / (v[2][1] - v[0][1])) * (v[2][2] - v[0][2])))]
        bottomTriangle(v[0], v[1], v3, f);
        topTriangle(v[1], v3, v[2], f);

def zbuf_draw(t, scale, offset, f):
    if len(t) != 3:
        print("NOT A TRIANGLE")
        exit(-1)

    # some scaling ...
    lt = []
    for c in t:
        lt.append( [ int(scale[0]*c[0]+offset[0]), int(scale[1]*c[1]+offset[1]), c[2] ] )

    drawTriangle(lt, f)

def zbuf_dump():
    # dump into image for visual test
    # convert -size 1024x768 -depth 8 gray:image.raw image.png  
    with open('image.raw', 'wb') as f:        
        for y in zbuf:
            for x in y:
                if x == None: f.write(bytes([0]))
                else:         f.write(bytes([255]))

import functools

def zbuf_eval():
    for line in zbuf:
        for pix in line:
            if pix != None:
                pix[1]["visible"] = True
        
# ============================== VRML/X3D parser ==============================

def parse_appearance(a):
    color = None
    
    #print("Appearance", a);
    for e in a:
        if e.tag == "Material":
            dc = e.attrib["diffuseColor"]
            c = dc.split()
            color = [ float(c[0]), float(c[1]), float(c[2]) ]
        else:
            print("Ignored:", a.tag);

    return color;
    
def parse_coordIndex(ci):
    c = [ ]
    coo = [ ]
    for l in ci.split():
        s = l.split()
        for ps in s:
            coo.append(int(ps))
            if len(coo) == 4:
                c.append(coo)
                coo = [ ]
    return c

def parse_point(pi):
    p = [ ]
    px = [ ]
    for l in pi.split():
        s = l.split()
        for ps in s:
            px.append(float(ps))
            if len(px) == 3:
                p.append(px)
                px = [ ]
    return p

def parse_coordinate(c):
    # print("Coordinate", c);
    pi = parse_point(c.attrib["point"]);
    # print("# points:", len(pi))
    for e in c:
        print("Ignored:", e.tag);

    return pi
        
def parse_indexedfaceset(i):
    # print("IndexedFaceSet", i);
    ci = parse_coordIndex(i.attrib["coordIndex"]);

    for e in i:
        if e.tag == "Coordinate":
            p = parse_coordinate(e)
        else:
            print("Ignored:", e.tag);

    # combine ci and p
    faceset = [ ]
    for c in ci:
        faceset.append( [ p[c[0]], p[c[1]], p[c[2]] ] )

    return faceset

def dotproduct(v1, v2):
    return sum((a*b) for a, b in zip(v1, v2))

def length(v):
    return math.sqrt(dotproduct(v, v))

def angle(v1, v2):
    a = dotproduct(v1, v2) / (length(v1) * length(v2))
    if a < -1: return math.pi;
    if a >  1: return 0;
    return math.acos(a)

# cross product of two vectors
def cross(a, b):
    c = [a[1]*b[2] - a[2]*b[1],
         a[2]*b[0] - a[0]*b[2],
         a[0]*b[1] - a[1]*b[0]]

    return c

def normal(f):
    a = [ f[1][0]-f[0][0], f[1][1]-f[0][1], f[1][2]-f[0][2] ]
    b = [ f[2][0]-f[0][0], f[2][1]-f[0][1], f[2][2]-f[0][2] ]
    return cross(a,b)
    
def parse_shape(s):
    shape = { }
    shape["name"] = s.attrib["DEF"]    
    # print("Shape", "#"+s.attrib["DEF"]);
    for e in s:
        if e.tag == "Appearance":
            shape["color"] = parse_appearance(e)
        elif e.tag == "IndexedFaceSet":
            shape["faces"] = parse_indexedfaceset(e)

            # we now have the complete faceset in fs
            # print("Faces:", len(shape["faces"]));
        else:
            print("Ignored:", e.tag);

    return shape
            
def parse_transform(t):
    print("Transform", str(t));

def face_distance(a):
    dist = 0
    for p in a["face"]: 
        dist = dist + p[2]
    return dist / len(a["face"])

def point_is_same(p0,p1):
    LIMIT = 0.0001
    if ((abs(p0[0]-p1[0]) < LIMIT) and
        (abs(p0[1]-p1[1]) < LIMIT) and
        (abs(p0[2]-p1[2]) < LIMIT)):
        return True

    return False

def fuse_polygons(af, bf):                
    # search for matching points
    for a in range(len(af["face"])):
        for b in range(len(bf["face"])):
            if point_is_same(af["face"][a], bf["face"][b]):
                an = 0;
                if a < len(af["face"])-1: an = a + 1
                bn = len(bf["face"])-1;
                if b > 0: bn = b - 1

                # check if "next" point also matches
                if point_is_same(af["face"][an], bf["face"][bn]):
                    # print("Same", a, an, "->", b, bn);
                    
                    # 0 ... a of first polygon, then
                    # b+1 ... bn-1 of second, then
                    # an ... amax of first

                    new_pts = [ ]
                    
                    # 0 ... a of first polygon
                    new_pts.extend(af["face"][0:a+1])
                    # b+1 ... bn-1 of second, this may wrap as b+1 may be > bn-1
                    if bn > b:
                        new_pts.extend(bf["face"][b+1:bn+1])
                    else:
                        new_pts.extend(bf["face"][b+1:])
                        new_pts.extend(bf["face"][:bn+1])
                    # an ... amax of first
                    if an > 0:
                        new_pts.extend(af["face"][an+1:])

                    if len(new_pts) < 4:
                        print("HUH???");
                        exit(1)

                    return new_pts

    return None

def fuse_faces(faces):
    fused = [ ]
    for f in faces:
        print(".", end=""); sys.stdout.flush()

        # compare current face with all faces already fused
        v = None
        r = None
        for fu in fused:
            if v == None:
                # check if face normals point into same direction
                if abs(angle(f["normal"], fu["normal"]))  < 0.0001:
                    r = fu
                    v = fuse_polygons(f, fu)

        # if this polygon could not be fused with any other then just
        # append it. Otherwise replace the polygon of the now fused
        # face
        if v == None:                
            fused.append(f)
        else:
            r["face"] = v

    print("")
    return fused
            
def parse_scene(s):
    scene = { }
    scene["shapes"] = [ ]
    
    # print("Scene", s);
    for e in s:
        if e.tag == "Shape":
            scene["shapes"].append(parse_shape(e))
        elif e.tag == "Transform":
            parse_transform(e)
        else:
            print("Ignored:", e.tag);

    # process scene
    print("# shapes: ", len(scene["shapes"]))

    # move all faces into one list
    faces = []
    id = 0
    for sh in scene["shapes"]:
        for fc in sh["faces"]:
            faces.append( { "id": id, "face": fc, "color": copy.copy(sh["color"]) } )
            id = id + 1

    # rotate
    r = [ 0.2,0.4,0 ]
    for f in faces:
        cn = [ ]
        for c in f["face"]:
            # rotate around x (horizontal)
            c1 = c[1] *  math.cos(r[0]) - c[2] * math.sin(r[0])
            c2 = c[1] *  math.sin(r[0]) + c[2] * math.cos(r[0])
            c = [ c[0], c1, c2 ]
            
            # rotate around y (vertical)
            c0 = c[0] *  math.cos(r[1]) + c[2] * math.sin(r[1])
            c2 = c[0] * -math.sin(r[1]) + c[2] * math.cos(r[1])
            c = [ c0, c[1], c2 ]
            
            # rotate around z (depth)
            c0 = c[0] *  math.cos(r[2]) - c[1] * math.sin(r[2])
            c1 = c[0] *  math.sin(r[2]) + c[1] * math.cos(r[2])
            c = [ c0, c1, c[2] ]
            
            cn.append(c)
        f["face"] = cn
            
    # create face normals
    for f in faces:
        f["normal"] = normal(f["face"])

    # remove all faces that are seen from the rear side
    oldf = faces
    faces = [ ]
    for f in oldf:
        # if f["normal"][2] > 0:
        if angle(f["normal"], [ 0,0,1]) < math.pi/2:
            faces.append(f)

    print("Faces facing away from user:", len(oldf) - len(faces), "visible:",  len(faces))
            
    # sort by depth
    faces = sorted(faces, key=face_distance)

    # perspective projection
    # TODO
    
    # determine min/max range
    min_max_range = [ [ 10000, -10000 ], [ 10000, -10000 ], [ 10000, -10000 ] ]
    for f in faces:
        for coo in f["face"]:
            if coo[0] < min_max_range[0][0]: min_max_range[0][0] = coo[0]
            if coo[0] > min_max_range[0][1]: min_max_range[0][1] = coo[0]
            if coo[1] < min_max_range[1][0]: min_max_range[1][0] = coo[1]
            if coo[1] > min_max_range[1][1]: min_max_range[1][1] = coo[1]
            if coo[2] < min_max_range[2][0]: min_max_range[2][0] = coo[2]
            if coo[2] > min_max_range[2][1]: min_max_range[2][1] = coo[2]

    print("Range X/Y/Z: ", min_max_range);

    scale = [ int(1000/(min_max_range[0][1] - min_max_range[0][0])), int(1000/(min_max_range[1][1] - min_max_range[1][0])) ]
    print("SCALE:", scale)
    offset = [ -min_max_range[0][0]*scale[0], -min_max_range[1][0]*scale[1] ]
    print("OFFSET:", offset)
    
    # draw into z buffer, at this point all faces are still
    # triangles
    for f in faces:
        f["visible"] = False
        zbuf_draw(f["face"], scale, offset, f)

    # dump into file to see if drawing works
    # zbuf_dump();

    # mark faces visible in final image
    zbuf_eval();

    # drop all invisble faces
    visible = []
    for f in faces:
        if f["visible"]:
            visible.append(f)

    print("Visible faces:", len(visible));
    faces = visible

    # depth sorting by zbuffer -> evaluate coverbuf
#    for i in list(coverbuf.keys()):
#        print("ID", i, "covers", coverbuf[i])

    # find all faces that don't cover anything
    newf = [ ]

    #    while len(faces) > 0:
    if False:
        print("Processing faces", len(faces))
        
        keepf = [ ]
   
        # do this as long as there are faces left
        for f in faces:
            # check if this face has an entry in the coverbuf
            # (the list containing which faces a face covers)
            if not f["id"] in coverbuf:
                print("does not cover anything", f["id"])
                newf.append(f)
                rem = [ ]
                for c in coverbuf:
                    if f["id"] in coverbuf[c]:
#                        print("C:", c, coverbuf[c])
                        coverbuf[c].remove(f["id"])
                        if len(coverbuf[c]) == 0:
                            rem.append(c)
#                        print("->", c, coverbuf[c])
                for c in rem:
#                    print("rmoving", c)
                    coverbuf.pop(c)
            else:
                keepf.append(f)

        faces = keepf
            
#    faces = newf
    
    # fuse faces
    print("Fusing");
    pfaces = len(faces)
    faces = fuse_faces(faces);
    print("\nFused:", len(faces))
    while len(faces) != pfaces:
        print("\nFusing");
        pfaces = len(faces)
        faces = fuse_faces(faces);
        print("\nFused:", len(faces))
        
    # shade with respect to viewing angle
    for f in faces:
        shade = 1.4*(1.0-(angle(f["normal"], [ 0,0,1]) / (math.pi/2)))
        # print("shade", shade)
        f["color"][0] = f["color"][0]*shade
        f["color"][1] = f["color"][1]*shade
        f["color"][2] = f["color"][2]*shade
        if f["color"][0] > 1: f["color"][0] = 1
        if f["color"][1] > 1: f["color"][1] = 1
        if f["color"][2] > 1: f["color"][2] = 1

    print("Drawing faces", len(faces))
        
    svg_document = pysvg.structure.Svg()
    shape_builder = pysvg.builders.ShapeBuilder()
    style_builder = pysvg.builders.StyleBuilder()
    style_builder.setStrokeLineJoin("round")
    
    # draw, svg size is 100x100
    offset = [ 80, 80 ]
    scale = 20.0
    for f in faces:
        pts = []
        # draw polygons
        for c in f["face"]:
            pts.append( (scale*c[0]+offset[0], -scale*c[1]+offset[1]) )

        ps = shape_builder.convertTupleArrayToPoints(pts)
        color = "#%02x%02x%02x" % (int(255*f["color"][0]), int(255*f["color"][1]), int(255*f["color"][2]) )
        poly = shape_builder.createPolygon(points=ps, strokewidth="0.1px", stroke = "black", fill = color)
        #        poly.set_style(style_builder.getStyle())
        svg_document.addElement(poly)
        # svg_document.addElement(shape_builder.createPolygon(points=ps, strokewidth="0.1px", stroke = color, fill = color))

        # stroke-linejoin:round
        
    svg_document.save("ftd.svg")
            
def parse_file(n):
    # try to load file, just bail out if anything goes wrong
    tree = et.parse(n)

    root = tree.getroot()
    if root.tag != "X3D":
        print("Expecting X3D")
        return False

    # parse_x3d(root)
    for e in root:
        if e.tag == "Scene":
            parse_scene(e)

    return True
    
def main():
    if len(sys.argv) != 2:
        print("No arguments given")
        return False
    
    print("Parsing", sys.argv[1], "...");
    if not parse_file(sys.argv[1]):
        return False

    return True

if __name__ == "__main__":
    if not main():
        exit(-1);
