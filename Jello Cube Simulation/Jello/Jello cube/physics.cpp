/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code


  Dhruv Rajvansh

*/

#include "jello.h"
#include "physics.h"




const double rLStruc = 1.0 / 7.0;
const double rLShearFace = sqrt(2.0) * rLStruc;
const double rLShearMainDiag = sqrt(3.0) * rLStruc;
const double rLBend = 2 * rLStruc;

double pLength(point L) {
    double le = sqrt((L).x * (L).x + (L).y * (L).y + (L).z * (L).z);
    if (le < 1e-7) return 0;
    return le;
}

double pDotProduct(point a, point b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

void pNormalize(point v)
{
    double length = pLength(v);
    v.x /= length;
    v.y /= length;
    v.z /= length;
}


point ComputHooks(struct point A, struct point B, double kh, double rl)
{
    point L ;
    pDIFFERENCE(A, B, L);
    point NormL;
    pCPY(L, NormL);
    pNormalize(NormL);
    double ScaleL;
    double lenL = pLength(L);
    double Scalae = (-kh) * (lenL - rl);
    point F;
    pMULTIPLY(NormL, Scalae, F);
    return F;
}

point ComputeDamping(struct point A, struct point B, struct point Va, struct point Vb, double kb) {

    point L = {0,0,0};
    pDIFFERENCE(A, B, L);
    point vel;
    pDIFFERENCE(Va, Vb, vel);
    point NormL;
    pCPY(L, NormL);
    pNormalize(NormL);
    double LengthyL = pLength(L);
    if (LengthyL < 1e-7) {
        point zero = { 0,0,0 };
        return zero;
    }
    double pr = -kb * pDotProduct(vel, L) / LengthyL;
    point F;
    pMULTIPLY(NormL, pr, F);
    return F;
}

point ComputeStructForce(struct world* jello, int i, int j, int k) {

    point temp = {0,0,0};

    if (i - 1 >= 0) {
        point hooks = ComputHooks(jello->p[i][j][k], jello->p[i - 1][j][k], jello->kElastic, rLStruc);
        point Damp = ComputeDamping(jello->p[i][j][k], jello->p[i - 1][j][k], jello->v[i][j][k], jello->v[i - 1][j][k], jello->dElastic);
        pSUM(temp, hooks, temp);
        pSUM(temp, Damp, temp);
    }
    if (i + 1 <= 7) {
        point hooks = ComputHooks(jello->p[i][j][k], jello->p[i + 1][j][k], jello->kElastic, rLStruc);
        point Damp = ComputeDamping(jello->p[i][j][k], jello->p[i + 1][j][k], jello->v[i][j][k], jello->v[i + 1][j][k], jello->dElastic);
        pSUM(temp, hooks, temp);
        pSUM(temp, Damp, temp);
    }
    if (j - 1 >= 0) {
        point hooks = ComputHooks(jello->p[i][j][k], jello->p[i][j - 1][k], jello->kElastic, rLStruc);
        point Damp = ComputeDamping(jello->p[i][j][k], jello->p[i][j - 1][k], jello->v[i][j][k], jello->v[i][j - 1][k], jello->dElastic);
        pSUM(temp, hooks, temp);
        pSUM(temp, Damp, temp);
    }
    if (j + 1 <= 7) {
        point hooks = ComputHooks(jello->p[i][j][k], jello->p[i][j + 1][k], jello->kElastic, rLStruc);
        point Damp = ComputeDamping(jello->p[i][j][k], jello->p[i][j + 1][k], jello->v[i][j][k], jello->v[i][j + 1][k], jello->dElastic);
        pSUM(temp, hooks, temp);
        pSUM(temp, Damp, temp);
    }
    if (k - 1 >= 0) {
        point hooks = ComputHooks(jello->p[i][j][k], jello->p[i][j][k - 1], jello->kElastic, rLStruc);
        point Damp = ComputeDamping(jello->p[i][j][k], jello->p[i][j][k - 1], jello->v[i][j][k], jello->v[i][j][k - 1], jello->dElastic);
        pSUM(temp, hooks, temp);
        pSUM(temp, Damp, temp);
    }
    if (k + 1 <= 7) {
        point hooks = ComputHooks(jello->p[i][j][k], jello->p[i][j][k + 1], jello->kElastic, rLStruc);
        point Damp = ComputeDamping(jello->p[i][j][k], jello->p[i][j][k + 1], jello->v[i][j][k], jello->v[i][j][k + 1], jello->dElastic);
        pSUM(temp, hooks, temp);
        pSUM(temp, Damp, temp);
    }
    return temp;
}

point ComputeBendForce(struct world *jello, int i, int j, int k) {
    point temp = {0,0,0};
    if (i - 2 >= 0) {
        point hooks = ComputHooks(jello->p[i][j][k], jello->p[i-2][j][k], jello->kElastic, rLBend);
        point Damp = ComputeDamping(jello->p[i][j][k], jello->p[i-2][j][k], jello->v[i][j][k], jello->v[i-2][j][k], jello->dElastic);

        pSUM(temp, hooks, temp);
        pMULTIPLY(Damp, 1.5, Damp);
        pSUM(temp, Damp, temp);
    }

    if (i + 2 <= 7) {
        point hooks = ComputHooks(jello->p[i][j][k], jello->p[i + 2][j][k], jello->kElastic, rLBend);
        point Damp = ComputeDamping(jello->p[i][j][k], jello->p[i + 2][j][k], jello->v[i][j][k], jello->v[i + 2][j][k], jello->dElastic);
        pSUM(temp, hooks, temp);
        pMULTIPLY(Damp, 1.5, Damp);
        pSUM(temp, Damp, temp);
    }
    if (j - 2 >= 0) {
        point hooks = ComputHooks(jello->p[i][j][k], jello->p[i][j-2][k], jello->kElastic, rLBend);
        point Damp = ComputeDamping(jello->p[i][j][k], jello->p[i][j-2][k], jello->v[i][j][k], jello->v[i][j-2][k], jello->dElastic);
        pSUM(temp, hooks, temp);
        pMULTIPLY(Damp, 1.5, Damp);
        pSUM(temp, Damp, temp);
    }

    if (j + 2 <= 7) {
        point hooks = ComputHooks(jello->p[i][j][k], jello->p[i][j+2][k], jello->kElastic, rLBend);
        point Damp = ComputeDamping(jello->p[i][j][k], jello->p[i][j+2][k], jello->v[i][j][k], jello->v[i][j+2][k], jello->dElastic);
        pSUM(temp, hooks, temp);
        pSUM(temp, Damp, temp);
    }
    if (k - 2 >= 0) {
        point hooks = ComputHooks(jello->p[i][j][k], jello->p[i][j][k-2], jello->kElastic, rLBend);
        point Damp = ComputeDamping(jello->p[i][j][k], jello->p[i][j][k-2], jello->v[i][j][k], jello->v[i][j][k-2], jello->dElastic);
        pSUM(temp, hooks, temp);
        pMULTIPLY(Damp, 1.5, Damp);
        pSUM(temp, Damp, temp);
    }

    if (k + 2 <= 7) {
        point hooks = ComputHooks(jello->p[i][j][k], jello->p[i][j][k+2], jello->kElastic, rLBend);
        point Damp = ComputeDamping(jello->p[i][j][k], jello->p[i][j][k+2], jello->v[i][j][k], jello->v[i][j][k+2], jello->dElastic);
        pSUM(temp, hooks, temp);
        pMULTIPLY(Damp, 1.5, Damp);
        pSUM(temp, Damp, temp);
    }
    return temp;
}

point ComputeShearForce(struct world* jello, int i, int j, int k) {
    point temp = {0,0,0};
    if ((j - 1 >= 0) && (k - 1 >= 0)) {
        point hooks = ComputHooks(jello->p[i][j][k], jello->p[i][j-1][k -1], jello->kElastic, rLShearFace);
        point Damp = ComputeDamping(jello->p[i][j][k], jello->p[i][j - 1][k - 1],jello->v[i][j][k], jello->v[i][j-1][k-1], jello->dElastic);
        pSUM(temp,hooks,temp);
        pMULTIPLY(Damp, 1.5, Damp);
        pSUM(temp, Damp, temp);
    }
    if ((j-1 >=0) && (k+1 <= 7)) {
        point hooks = ComputHooks(jello->p[i][j][k], jello->p[i][j - 1][k + 1], jello->kElastic, rLShearFace);
        point Damp = ComputeDamping(jello->p[i][j][k], jello->p[i][j - 1][k + 1], jello->v[i][j][k], jello->v[i][j - 1][k + 1], jello->dElastic);
        pSUM(temp, hooks, temp);
        pMULTIPLY(Damp, 1.5, Damp);
        pSUM(temp, Damp, temp);
    }
    if ((j+1 <= 7)&&(k-1 >=0)) {
        point hooks = ComputHooks(jello->p[i][j][k], jello->p[i][j + 1][k - 1], jello->kElastic, rLShearFace);
        point Damp = ComputeDamping(jello->p[i][j][k], jello->p[i][j + 1][k - 1], jello->v[i][j][k], jello->v[i][j + 1][k - 1], jello->dElastic);
        pSUM(temp, hooks, temp);
        pMULTIPLY(Damp, 1.5, Damp);
        pSUM(temp, Damp, temp);
    }
    if ((j+1  <=7) && (k+1 <= 7)) {
        point hooks = ComputHooks(jello->p[i][j][k], jello->p[i][j + 1][k + 1], jello->kElastic, rLShearFace);
        point Damp = ComputeDamping(jello->p[i][j][k], jello->p[i][j + 1][k + 1], jello->v[i][j][k], jello->v[i][j + 1][k + 1], jello->dElastic);
        pSUM(temp, hooks, temp);
        pMULTIPLY(Damp, 1.5, Damp);
        pSUM(temp, Damp, temp);
    }
    if ((i - 1 >= 0) && (k - 1 >= 0)) {
        point hooks = ComputHooks(jello->p[i][j][k], jello->p[i-1][j][k - 1], jello->kElastic, rLShearFace);
        point Damp = ComputeDamping(jello->p[i][j][k], jello->p[i-1][j][k - 1], jello->v[i][j][k], jello->v[i-1][j][k - 1], jello->dElastic);
        pSUM(temp, hooks, temp);
        pMULTIPLY(Damp, 1.5, Damp);
        pSUM(temp, Damp, temp);
    }
    if ((i - 1 >= 0) && (k + 1 <= 7)) {
        point hooks = ComputHooks(jello->p[i][j][k], jello->p[i-1][j][k + 1], jello->kElastic, rLShearFace);
        point Damp = ComputeDamping(jello->p[i][j][k], jello->p[i-1][j][k + 1], jello->v[i][j][k], jello->v[i-1][j][k + 1], jello->dElastic);
        pSUM(temp, hooks, temp);
        pMULTIPLY(Damp, 1.5, Damp);
        pSUM(temp, Damp, temp);
    }
    if ((i + 1 <= 7) && (k - 1 >= 0)) {
        point hooks = ComputHooks(jello->p[i][j][k], jello->p[i + 1][j][k - 1], jello->kElastic, rLShearFace);
        point Damp = ComputeDamping(jello->p[i][j][k], jello->p[i + 1][j][k - 1], jello->v[i][j][k], jello->v[i + 1][j][k - 1], jello->dElastic);
        pSUM(temp, hooks, temp);
        pMULTIPLY(Damp, 1.5, Damp);
        pSUM(temp, Damp, temp);
    }
    if ((i + 1 <= 7) && (k + 1 <= 7)) {
        point hooks = ComputHooks(jello->p[i][j][k], jello->p[i + 1][j][k + 1], jello->kElastic, rLShearFace);
        point Damp = ComputeDamping(jello->p[i][j][k], jello->p[i + 1][j][k + 1], jello->v[i][j][k], jello->v[i + 1][j][k + 1], jello->dElastic);
        pSUM(temp, hooks, temp);
        pMULTIPLY(Damp, 1.5, Damp);
        pSUM(temp, Damp, temp);
    }
    if ((i - 1 >= 0) && (j - 1 >= 0)) // If A has left, bottom face diagonal neighbor
    {
        point hooks =  ComputHooks(jello->p[i][j][k], jello->p[i - 1][j - 1][k], jello->kElastic, rLShearFace);
        pSUM(temp, hooks, temp);
        point Damp = ComputeDamping(jello->p[i][j][k], jello->p[i - 1][j - 1][k], jello->v[i][j][k], jello->v[i - 1][j - 1][k], jello->dElastic);
        pMULTIPLY(Damp, 1.5, Damp);
        pSUM(Damp, temp, temp);
        if (k - 1 >= 0)  // For left, bottom, back 3D cube main diagonal
        {
            point hooks = ComputHooks(jello->p[i][j][k], jello->p[i - 1][j - 1][k - 1], jello->kElastic, rLShearMainDiag);
            pSUM(hooks, temp, temp);
            point Damp = ComputeDamping(jello->p[i][j][k], jello->p[i - 1][j - 1][k - 1], jello->v[i][j][k], jello->v[i - 1][j - 1][k - 1], jello->dElastic);
            pMULTIPLY(Damp, 1.5, Damp);
            pSUM(temp, Damp, temp);
        }
        if (k + 1 <= 7)  // For left, bottom, front 3D cube main diagonal
        {
            point hooks = ComputHooks(jello->p[i][j][k], jello->p[i - 1][j - 1][k + 1], jello->kElastic, rLShearMainDiag);
            pSUM(hooks, temp, temp);
            point Damp = ComputeDamping(jello->p[i][j][k], jello->p[i - 1][j - 1][k + 1], jello->v[i][j][k], jello->v[i - 1][j - 1][k + 1], jello->dElastic);
            pMULTIPLY(Damp, 1.5, Damp);
            pSUM(Damp, temp, temp);
        }
    }
    if ((i - 1 >= 0) && (j + 1 <= 7)) // If A has left, top face diagonal neighbor
    {
        point hooks = ComputHooks(jello->p[i][j][k], jello->p[i - 1][j + 1][k], jello->kElastic, rLShearFace);
        pSUM(temp, hooks, temp);
        point Damp = ComputeDamping(jello->p[i][j][k], jello->p[i - 1][j + 1][k], jello->v[i][j][k], jello->v[i - 1][j + 1][k], jello->dElastic);
        
        pSUM(Damp, temp, temp);
        if (k - 1 >= 0)  // For left, top, back 3D cube main diagonal
        {
            point hooks = ComputHooks(jello->p[i][j][k], jello->p[i - 1][j + 1][k - 1], jello->kElastic, rLShearMainDiag);
            pSUM(hooks, temp, temp);
            point Damp = ComputeDamping(jello->p[i][j][k], jello->p[i - 1][j + 1][k - 1], jello->v[i][j][k], jello->v[i - 1][j + 1][k - 1], jello->dElastic);
            pMULTIPLY(Damp, 1.5, Damp);
            pSUM(Damp, temp, temp);
        }
        if (k + 1 <= 7)  // For left, top, front 3D cube main diagonal
        {
            point hooks = ComputHooks(jello->p[i][j][k], jello->p[i - 1][j + 1][k + 1], jello->kElastic, rLShearMainDiag);
            pSUM(hooks, temp, temp);
            point Damp = ComputeDamping(jello->p[i][j][k], jello->p[i - 1][j + 1][k + 1], jello->v[i][j][k], jello->v[i - 1][j + 1][k + 1], jello->dElastic);
            pMULTIPLY(Damp, 1.5, Damp);
            pSUM(Damp, temp, temp);
        }
    }
    if ((i + 1 <= 7) && (j - 1 >= 0)) // If A has right, bottom face diagonal neighbor
    {
        point hooks= ComputHooks(jello->p[i][j][k], jello->p[i + 1][j - 1][k], jello->kElastic, rLShearFace);
        pSUM(hooks, temp, temp);
        point Damp = ComputeDamping(jello->p[i][j][k], jello->p[i + 1][j - 1][k], jello->v[i][j][k], jello->v[i + 1][j - 1][k], jello->dElastic);
        
        pSUM(Damp, temp, temp);
        if (k - 1 >= 0)  // For right, bottom, back 3D cube main diagonal
        {
            point hooks = ComputHooks(jello->p[i][j][k], jello->p[i + 1][j - 1][k - 1], jello->kElastic, rLShearMainDiag);
            pSUM(hooks, temp, temp);
            point Damp = ComputeDamping(jello->p[i][j][k], jello->p[i + 1][j - 1][k - 1], jello->v[i][j][k], jello->v[i + 1][j - 1][k - 1], jello->dElastic);
            pMULTIPLY(Damp, 1.5, Damp);
            pSUM(Damp, temp, temp);
        }
        if (k + 1 <= 7)  // For right, bottom, front 3D cube main diagonal
        {
            point hooks = ComputHooks(jello->p[i][j][k], jello->p[i + 1][j - 1][k + 1], jello->kElastic, rLShearMainDiag);
            pSUM(hooks, temp, temp);
            point Damp = ComputeDamping(jello->p[i][j][k], jello->p[i + 1][j - 1][k + 1], jello->v[i][j][k], jello->v[i + 1][j - 1][k + 1], jello->dElastic);
            pMULTIPLY(Damp, 1.5, Damp);
            pSUM(Damp, temp, temp);
        }
    }
    if ((i + 1 <= 7) && (j + 1 <= 7)) // If A has right, top face diagonal neighbor
    {
        point hooks = ComputHooks(jello->p[i][j][k], jello->p[i + 1][j + 1][k], jello->kElastic, rLShearFace);
        pSUM(hooks, temp, temp);
        point Damp = ComputeDamping(jello->p[i][j][k], jello->p[i + 1][j + 1][k], jello->v[i][j][k], jello->v[i + 1][j + 1][k], jello->dElastic);
        
        pSUM(Damp, temp, temp);
        if (k - 1 >= 0)  // For right, top, back 3D cube main diagonal
        {
            point hooks = ComputHooks(jello->p[i][j][k], jello->p[i + 1][j + 1][k - 1], jello->kElastic, rLShearMainDiag);
            pSUM(hooks, temp, temp);
            point Damp = ComputeDamping(jello->p[i][j][k], jello->p[i + 1][j + 1][k - 1], jello->v[i][j][k], jello->v[i + 1][j + 1][k - 1], jello->dElastic);
            pMULTIPLY(Damp, 1.5, Damp);
            pSUM(Damp, temp, temp);
        }
        if (k + 1 <= 7)  // For right, top, front 3D cube main diagonal
        {
            point hooks = ComputHooks(jello->p[i][j][k], jello->p[i + 1][j + 1][k + 1], jello->kElastic, rLShearMainDiag);
            pSUM(hooks, temp, temp);
            point Damp = ComputeDamping(jello->p[i][j][k], jello->p[i + 1][j + 1][k + 1], jello->v[i][j][k], jello->v[i + 1][j + 1][k + 1], jello->dElastic);
            pMULTIPLY(Damp, 1.5,Damp);
            pSUM(Damp, temp, temp);
        }
    }
    return temp;
}

// Computes collision spring forces on mass point A from collision with bounding box wall using Penalty method
point computeCollisionForce(struct world* jello, int i, int j, int k)
{
    point temp = {0,0,0};
    point vObstacle = {0,0,0};
 
    point cp[6];
    double mag; // Magnitude will be proportional to the amount of penetration
    pMAKE(2, jello->p[i][j][k].y, jello->p[i][j][k].z, cp[0]); // initializing contact points for each wall
    pMAKE(-2, jello->p[i][j][k].y, jello->p[i][j][k].z, cp[1]);
    pMAKE(jello->p[i][j][k].x, 2, jello->p[i][j][k].z, cp[2]);
    pMAKE(jello->p[i][j][k].x, -2, jello->p[i][j][k].z, cp[3]);
    pMAKE(jello->p[i][j][k].x, jello->p[i][j][k].y, 2, cp[4]);
    pMAKE(jello->p[i][j][k].x, jello->p[i][j][k].y, -2, cp[5]);
    if (jello->p[i][j][k].x > 2)
    {
        mag = fmax(fabs(jello->p[i][j][k].x) - 2, 0);
        point hooks = ComputHooks(jello->p[i][j][k], cp[0], (jello->kCollision ), 0);
        pSUM(hooks, temp, temp);
        point Damp = ComputeDamping(jello->p[i][j][k], cp[0], jello->v[i][j][k], vObstacle, (jello->dCollision));
        pMULTIPLY(Damp, 3, Damp);
        pSUM(Damp, temp, temp);
    }
    if (jello->p[i][j][k].x < -2)
    {
        mag = fmax(fabs(jello->p[i][j][k].x) - 2, 0);
        point hooks = ComputHooks(jello->p[i][j][k], cp[1], (jello->kCollision), 0);
        pSUM(hooks, temp, temp);
        point Damp = ComputeDamping(jello->p[i][j][k], cp[1], jello->v[i][j][k], vObstacle, (jello->dCollision ));
        pMULTIPLY(Damp, 3, Damp);
        pSUM(Damp, temp, temp);
    }
    if (jello->p[i][j][k].y > 2)
    {
        mag = fmax(fabs(jello->p[i][j][k].y) - 2, 0);
        point hooks = ComputHooks(jello->p[i][j][k], cp[2], (jello->kCollision ), 0);
        pSUM(hooks, temp, temp);
        point Damp = ComputeDamping(jello->p[i][j][k], cp[2], jello->v[i][j][k], vObstacle, (jello->dCollision));
        pMULTIPLY(Damp, 3, Damp);
        pSUM(Damp, temp, temp);
    }
    if (jello->p[i][j][k].y < -2)
    {
        mag = fmax(fabs(jello->p[i][j][k].y) - 2, 0);
        point hooks = ComputHooks(jello->p[i][j][k], cp[3], (jello->kCollision), 0);
        pSUM(hooks, temp, temp);
        point Damp = ComputeDamping(jello->p[i][j][k], cp[3], jello->v[i][j][k], vObstacle, (jello->dCollision ));
        pMULTIPLY(Damp, 3, Damp);
        pSUM(Damp, temp, temp);
    }
    if (jello->p[i][j][k].z > 2)
    {
        mag = fmax(fabs(jello->p[i][j][k].z) - 2, 0);
        point hooks = ComputHooks(jello->p[i][j][k], cp[4], (jello->kCollision * mag), 0);
        pSUM(hooks, temp, temp);
        point Damp = ComputeDamping(jello->p[i][j][k], cp[4], jello->v[i][j][k], vObstacle, (jello->dCollision));
        pMULTIPLY(Damp, 3, Damp);
        pSUM(Damp, temp, temp);
    }
    if (jello->p[i][j][k].z < -2)
    {
        mag = fmax(fabs(jello->p[i][j][k].z) - 2, 0);
        point hooks = ComputHooks(jello->p[i][j][k], cp[5], (jello->kCollision), 0);
        pSUM(hooks, temp, temp);
        point Damp = ComputeDamping(jello->p[i][j][k], cp[5], jello->v[i][j][k], vObstacle, (jello->dCollision));
        pMULTIPLY(Damp, 3, Damp);
        pSUM(Damp, temp, temp);
    }
    return temp;
}

point computeExternalForce(struct world* jello, int i, int j, int k)
{

    point temp = {0,0,0};
    double ffVoxelLen = 4.0 / (jello->resolution - 1);  // Length of a voxel in the force field
 
    int X = int((jello->p[i][j][k].x + 2) / ffVoxelLen);
    int Y = int((jello->p[i][j][k].y + 2) / ffVoxelLen);
    int Z = int((jello->p[i][j][k].z + 2) / ffVoxelLen);
    // Checking if (X,Y,Z) is within bounding box (0,0,0) to (resolution-1,resolution-1,resolution-1) in force field coordinate system
    if (X >= (jello->resolution - 1))
        X = (jello->resolution -2) ;
    else if (X < 0)
        X = 0;
    if (Y >= (jello->resolution - 1))
        Y = (jello->resolution - 2);
    else if (Y < 0)
        Y = 0;
    if (Z >= (jello->resolution - 1))
        Z = (jello->resolution - 2) ;
    else if (Z < 0)
        Z = 0;

    // Forces at 8 corners of voxel with mass point A
    point F000, F001, F010, F011, F100, F101, F110, F111;
    pCPY(jello->forceField[X * jello->resolution * jello->resolution + Y * jello->resolution + Z], F000);
    pCPY(jello->forceField[X * jello->resolution * jello->resolution + Y * jello->resolution + (Z + 1)], F001);
    pCPY(jello->forceField[X * jello->resolution * jello->resolution + (Y + 1) * jello->resolution + Z], F010);
    pCPY(jello->forceField[X * jello->resolution * jello->resolution + (Y + 1) * jello->resolution + (Z + 1)], F011);
    pCPY(jello->forceField[(X + 1) * jello->resolution * jello->resolution + Y * jello->resolution + Z], F100);
    pCPY(jello->forceField[(X + 1) * jello->resolution * jello->resolution + Y * jello->resolution + (Z + 1)], F101);
    pCPY(jello->forceField[(X + 1) * jello->resolution * jello->resolution + (Y + 1) * jello->resolution + Z], F110);
    pCPY(jello->forceField[(X + 1) * jello->resolution * jello->resolution + (Y + 1) * jello->resolution + (Z + 1)], F111);

    // Trilinear interpolation - 
    // Converting force field (X,Y,Z) position of mass point A back to world space and then computing [(i,j,k) - worldspaceA]/ffVoxelLen
    point A;
    A.x = (jello->p[i][j][k].x - (X * ffVoxelLen - 2)) / ffVoxelLen;
    A.y = (jello->p[i][j][k].y - (Y * ffVoxelLen - 2)) / ffVoxelLen;
    A.z = (jello->p[i][j][k].z - (Z * ffVoxelLen - 2)) / ffVoxelLen;
    // Computing the contribution of force at voxel corners to mass point A using Convex Combination
    pMULTIPLY(F000, ((1 - A.x) * (1 - A.y) * (1 - A.z)), F000);
    pSUM(F000, temp, temp);
    pMULTIPLY(F001, ((1 - A.x) * (1 - A.y) * A.z), F001);
    pSUM(F001, temp, temp);
    pMULTIPLY(F010, ((1 - A.x) * A.y * (1 - A.z)), F010);
    pSUM(F010, temp, temp);
    pMULTIPLY(F011, ((1 - A.x) * A.y * A.z), F011);
    pSUM(F011, temp, temp);
    pMULTIPLY(F100, (A.x * (1 - A.y) * (1 - A.z)), F100);
    pSUM(F100, temp, temp);
    pMULTIPLY(F101, (A.x * (1 - A.y) * A.z), F101);
    pSUM(F101, temp, temp);
    pMULTIPLY(F110, (A.x * A.y * (1 - A.z)), F110);
    pSUM(F110, temp,temp);
    pMULTIPLY(F111, (A.x * A.y * A.z), F111);
    pSUM(F111, temp, temp);

    return temp;
}

bool CheckCollistionToPlane(point p, point n, double d)
{
    double distanceToPlane = (n.x * p.x) + (n.y * p.y) + (n.z * p.z) + d;
    return distanceToPlane <= 0;
}

point CollideWithPlane(struct world* jello, int i, int j, int k) {
    point temp = { 0,0,0 };
    point N;
    pMAKE(jello->a, jello->b, jello->c, N);
    pNormalize(N);  // Normalize the plane normal
    double distance = jello->d;

    if (CheckCollistionToPlane(jello->p[i][j][k], N, distance)) {
        double penetrationDepth = -(N.x * jello->p[i][j][k].x + N.y * jello->p[i][j][k].y + N.z * jello->p[i][j][k].z + distance);

        point collisionPoint;
        collisionPoint.x = jello->p[i][j][k].x + N.x * penetrationDepth;
        collisionPoint.y = jello->p[i][j][k].y + N.y * penetrationDepth;
        collisionPoint.z = jello->p[i][j][k].z + N.z * penetrationDepth;

        point hooks = ComputHooks(jello->p[i][j][k], collisionPoint, jello->kCollision, 0);
        point Damp = ComputeDamping(jello->p[i][j][k], collisionPoint, jello->v[i][j][k],  {0, 0, 0 }, jello->dCollision);

        pMULTIPLY(Damp, 3, Damp);
        pSUM(hooks, temp, temp);
        pSUM(Damp, temp, temp);
    }

    return temp;
}


void computeAcceleration(struct world * jello, struct point a[8][8][8])
{
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            for (int k = 0; k < 8; k++) {
                point temp = { 0,0,0 };

                // Accumulate all forces
                point structForce = ComputeStructForce(jello, i, j, k);
                pSUM(temp, structForce, temp);

                point shearForce = ComputeShearForce(jello, i, j, k);
                pSUM(temp, shearForce, temp);

                point bendForce = ComputeBendForce(jello, i, j, k);
                pSUM(temp, bendForce, temp);

                point collisionForce = computeCollisionForce(jello, i, j, k);
                pSUM(temp, collisionForce, temp);

                point externalForce = computeExternalForce(jello, i, j, k);
                pSUM(temp, externalForce, temp);

                if (jello->incPlanePresent == 1) {
                    point planeforce = CollideWithPlane(jello,i,j,k);
                    pSUM(temp,planeforce,temp);
                }
                
                // Convert total force to acceleration (F = ma)
                pMULTIPLY(temp, 1.0f / jello->mass, a[i][j][k]);
            }
        }
    }
}

/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world * jello)
{
  int i,j,k;
  point a[8][8][8];

  computeAcceleration(jello, a);
  
  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
        jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
        jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
        jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;
        jello->v[i][j][k].x += jello->dt * a[i][j][k].x;
        jello->v[i][j][k].y += jello->dt * a[i][j][k].y;
        jello->v[i][j][k].z += jello->dt * a[i][j][k].z;

      }
}

/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */
void RK4(struct world * jello)
{
  point F1p[8][8][8], F1v[8][8][8], 
        F2p[8][8][8], F2v[8][8][8],
        F3p[8][8][8], F3v[8][8][8],
        F4p[8][8][8], F4v[8][8][8];

  point a[8][8][8];


  struct world buffer;

  int i,j,k;

  buffer = *jello; // make a copy of jello

  computeAcceleration(jello, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         pMULTIPLY(jello->v[i][j][k],jello->dt,F1p[i][j][k]);
         pMULTIPLY(a[i][j][k],jello->dt,F1v[i][j][k]);
         pMULTIPLY(F1p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F1v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F2p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F2p[i][j][k]);
         // F2v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F2v[i][j][k]);
         pMULTIPLY(F2p[i][j][k],0.5,buffer.p[i][j][k]);
         pMULTIPLY(F2v[i][j][k],0.5,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }

  computeAcceleration(&buffer, a);

  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F3p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F3v[i][j][k]);
         pMULTIPLY(F3p[i][j][k],1.0,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],1.0,buffer.v[i][j][k]);
         pSUM(jello->p[i][j][k],buffer.p[i][j][k],buffer.p[i][j][k]);
         pSUM(jello->v[i][j][k],buffer.v[i][j][k],buffer.v[i][j][k]);
      }
         
  computeAcceleration(&buffer, a);


  for (i=0; i<=7; i++)
    for (j=0; j<=7; j++)
      for (k=0; k<=7; k++)
      {
         // F3p = dt * buffer.v;
         pMULTIPLY(buffer.v[i][j][k],jello->dt,F4p[i][j][k]);
         // F3v = dt * a(buffer.p,buffer.v);     
         pMULTIPLY(a[i][j][k],jello->dt,F4v[i][j][k]);

         pMULTIPLY(F2p[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3p[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1p[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4p[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->p[i][j][k],jello->p[i][j][k]);

         pMULTIPLY(F2v[i][j][k],2,buffer.p[i][j][k]);
         pMULTIPLY(F3v[i][j][k],2,buffer.v[i][j][k]);
         pSUM(buffer.p[i][j][k],buffer.v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F1v[i][j][k],buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],F4v[i][j][k],buffer.p[i][j][k]);
         pMULTIPLY(buffer.p[i][j][k],1.0 / 6,buffer.p[i][j][k]);
         pSUM(buffer.p[i][j][k],jello->v[i][j][k],jello->v[i][j][k]);
      }

  return;  
}
