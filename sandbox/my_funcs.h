/*
static inline
coord embed_gradient (Point point, vector u, coord b, coord n)
{
  coord dudn;
  foreach_dimension() {
    bool dirichlet;
    double vb = u.x.boundary[embed] (point, point, u.x, &dirichlet);
    if (dirichlet) {
      double val;
      dudn.x = dirichlet_gradient (point, u.x, cs, n, b, vb, &val);
      dudn.x += u.x[]*val; // For pathological situations
    }
    else // Neumann
      dudn.x = vb;
    if (dudn.x == nodata)
      dudn.x = 0.;
  }
  return dudn;
}

double embed_interpolate_3d (Point point, scalar s, coord b)
{
  int i = sign(p.x), j = sign(p.y);
#if dimension == 2
  if (cs[i] && cs[0,j] && cs[i,j] &&
      (emerged || (csm1[i] && csm1[0,j] && csm1[i,j])))
    // bilinear interpolation when all neighbors are defined
    return ((s[]*(1. - fabs(p.x)) + s[i]*fabs(p.x))*(1. - fabs(p.y)) + 
	    (s[0,j]*(1. - fabs(p.x)) + s[i,j]*fabs(p.x))*fabs(p.y));
#else // dimension == 3
  int k = sign(p.z);
  if (cs[i,0,0] && cs[0,j,0] && cs[i,j,0] &&
      cs[0,0,k] && cs[i,0,k] && cs[0,j,k] && cs[i,j,k] &&
      (emerged || (csm1[i,0,0] && csm1[0,j,0] && csm1[i,j,0] &&
		   csm1[0,0,k] && csm1[i,0,k] && csm1[0,j,k] && csm1[i,j,k]))) {
    double val_0, val_k;
    // bilinear interpolation in x-y-planes when all neighbors are defined
    val_0 = (s[0,0,0]*(1. - fabs(p.x)) + s[i,0,0]*fabs(p.x))*(1. - fabs(p.y)) +
      (s[0,j,0]*(1. - fabs(p.x)) + s[i,j,0]*fabs(p.x))*fabs(p.y);
    val_k = (s[0,0,k]*(1. - fabs(p.x)) + s[i,0,k]*fabs(p.x))*(1. - fabs(p.y)) +
      (s[0,j,k]*(1. - fabs(p.x)) + s[i,j,k]*fabs(p.x))*fabs(p.y);
    // trilinear interpolation when all neighbors are defined
    return (val_0*(1. - fabs(p.z)) + val_k*fabs(p.z));
  }
#endif 
  else {
    // linear interpolation with gradients biased toward the
    // cells which are defined
    double val = s[];
    foreach_dimension() {
      int i = sign(p.x);
      if (cs[i] &&
	  (emerged || (csm1[] && csm1[i])))
	val += fabs(p.x)*(s[i] - s[]);
      else if (cs[-i] &&
	       (emerged || (csm1[] && csm1[-i])))
	val += fabs(p.x)*(s[] - s[-i]);
    }
    return val;
  }
}
*/

double embed_interpolate_3d_my (Point point, scalar s, coord p)
{
  int i = sign(p.x), j = sign(p.y);
#if dimension == 2
  if (cs[i] && cs[0,j] && cs[i,j]) {
    // bilinear interpolation when all neighbors are defined
    return ((s[]*(1. - fabs(p.x)) + s[i]*fabs(p.x))*(1. - fabs(p.y)) + 
	    (s[0,j]*(1. - fabs(p.x)) + s[i,j]*fabs(p.x))*fabs(p.y));
  }  
#else // dimension == 3
  int k = sign(p.z);
  if (cs[i,0,0] && cs[0,j,0] && cs[i,j,0] &&
      cs[0,0,k] && cs[i,0,k] && cs[0,j,k] && cs[i,j,k]) {
    double val_0, val_k;
    // bilinear interpolation in x-y-planes when all neighbors are defined
    val_0 = (s[0,0,0]*(1. - fabs(p.x)) + s[i,0,0]*fabs(p.x))*(1. - fabs(p.y)) +
      (s[0,j,0]*(1. - fabs(p.x)) + s[i,j,0]*fabs(p.x))*fabs(p.y);
    val_k = (s[0,0,k]*(1. - fabs(p.x)) + s[i,0,k]*fabs(p.x))*(1. - fabs(p.y)) +
      (s[0,j,k]*(1. - fabs(p.x)) + s[i,j,k]*fabs(p.x))*fabs(p.y);
    // trilinear interpolation when all neighbors are defined
    return (val_0*(1. - fabs(p.z)) + val_k*fabs(p.z));
  }
#endif 
  else {
    // linear interpolation with gradients biased toward the
    // cells which are defined
    double val = s[];
    foreach_dimension() {
      int i = sign(p.x);
      if (cs[i])
	val += fabs(p.x)*(s[i] - s[]);
      else if (cs[-i])
	val += fabs(p.x)*(s[] - s[-i]);
    }
    return val;
  }
}

trace
void embed_force_3d (scalar p, vector u, face vector mu, coord * Fp, coord * Fmu)
{
  coord Fps = {0}, Fmus = {0};
  foreach (reduction(+:Fps) reduction(+:Fmus)) {
    if (cs[] > 0. && cs[] < 1.) {
      coord n, b;
      double area = embed_geometry (point, &b, &n);
      area *= pow (Delta, dimension - 1);
      //double Fn = area*embed_interpolate (point, p, b);
      double Fn = area*embed_interpolate_3d_my (point, p, b);
      foreach_dimension() {
	Fps.x += Fn*n.x;
      }
      if (constant(mu.x) != 0.) {
	double mua = 0., fa = 0.;
	foreach_dimension() {
	  mua += mu.x[] + mu.x[1];
	  fa  += fs.x[] + fs.x[1];
	}
	mua /= (fa + SEPS);

	coord dudn = embed_gradient (point, u, b, n);
#if dimension == 2
	foreach_dimension()
	  Fmus.x -= area*mua*(dudn.x*(sq (n.x) + 1.) +
			      dudn.y*n.x*n.y);
#else // dimension == 3
	foreach_dimension()
	  Fmus.x -= area*mua*(dudn.x*(sq (n.x) + 1.) +
			      dudn.y*n.x*n.y +
			      dudn.z*n.x*n.z);
#endif // dimension
      }
    }
  }
  *Fp = Fps; *Fmu = Fmus;
}
