# Jello

Soft-body simulation extension for the renderer.

### Author: Matěj Fencl
### Category: Simulation / SolidObject
### Namespace: MatejFencl_Jello
### Class
Jello  
JelloManager  
JelloMesh : DefaultSceneNode, ISolid, ITimeDependent
### File: Jello.cs
### ITimeDependent: Yes

## Usage
### Preprocessing
Set floor collision height

    JelloManager.SetFloor(Y);

Create Jello using

    Jello jello = JelloManager.RegisterJello(ID, XS, YS, ZS, Scale, XOrigin, YOrigin, ZOrigin, Softness);

`ID` is unique identifier for this Jello object, used for JelloMesh.  
This will create Jello block consisting of `XS` * `YS` * `ZS` cube segments sized `Scale` (Recomended size is max cca. 16 cubic units due to perfomance and memory usage of BezierSurface)  
`XOrigin`, `YOrigin`, `ZOrigin` is center of the block.  
`Softness` is how soft/wobbly will the jello be. Default is 1, going too low will result in unstable simulation and setting the value too high will result in unusable blob of further unidentifiable shape (somewhere between 1 and 10 is ok).  

Add default forces for the Jello object

    jello.AddForce(XOrigin, YOrigin, ZOrigin, XForce, YForce, ZForce, Radius);
and then simulate the object

    jello.Simulate(TimeStart, TimeEnd, FPS);
### Rendering
Create new JelloMesh object

    JelloMesh jelloMesh = new JelloMesh(ID);
Use `ID` of Jello object created using preprocess stage.  
JelloMesh is regular Solid object and therefore can be added to scene as any other shape and will be animated.