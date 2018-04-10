// Animation time interval:
if ( outParam != null )
{
  outParam[ "Start" ] = 0.0;
  outParam[ "End" ]   = 20.0;
}

// CSG scene:
CSGInnerNode root = new CSGInnerNode( SetOperation.Union );
root.SetAttribute( PropertyName.REFLECTANCE_MODEL, new PhongModel() );
root.SetAttribute( PropertyName.MATERIAL, new PhongMaterial( new double[] { 1.0, 0.8, 0.1 }, 0.1, 0.6, 0.4, 128 ) );
scene.Intersectable = root;

// Background color:
scene.BackgroundColor = new double[] { 0.0, 0.05, 0.07 };

// Camera:
AnimatedCamera cam = new AnimatedCamera( new Vector3d( 0.7, -0.4,  0.0 ),
                                         new Vector3d( 0.7,  0.8, -6.0 ),
                                         50.0 );
cam.End = 20.0; // one complete turn takes 20.0 seconds
AnimatedRayScene ascene = scene as AnimatedRayScene;
if ( ascene != null )
  ascene.End = 20.0;
scene.Camera  = cam;

// Light sources:
scene.Sources = new LinkedList<ILightSource>();
scene.Sources.Add( new AmbientLightSource( 0.8 ) );
scene.Sources.Add( new PointLightSource( new Vector3d( -5.0, 4.0, -3.0 ), 1.2 ) );

// --- NODE DEFINITIONS ----------------------------------------------------

// Transparent sphere:
Sphere s;
s = new Sphere();
PhongMaterial pm = new PhongMaterial( new double[] { 0.0, 0.2, 0.1 }, 0.03, 0.03, 0.08, 128 );
pm.n  = 1.6;
pm.Kt = 0.9;
s.SetAttribute( PropertyName.MATERIAL, pm );
root.InsertChild( s, Matrix4d.Identity );

// Opaque sphere:
s = new Sphere();
root.InsertChild( s, Matrix4d.Scale( 1.2 ) * Matrix4d.CreateTranslation( 1.5, 0.2, 2.4 ) );

// Infinite plane with checker:
Plane pl = new Plane();
pl.SetAttribute( PropertyName.COLOR, new double[] { 0.0, 0.15, 0.0 } );
pl.SetAttribute( PropertyName.TEXTURE, new CheckerTexture( 0.6, 0.6, new double[] { 1.0, 1.0, 1.0 } ) );
root.InsertChild( pl, Matrix4d.RotateX( -MathHelper.PiOver2 ) * Matrix4d.CreateTranslation( 0.0, -1.0, 0.0 ) );
