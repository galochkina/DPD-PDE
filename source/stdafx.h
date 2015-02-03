#ifndef WINVER              // Allow use of features specific to Windows XP or later.
//#define WINVER 0x0501       // Change this to the appropriate value to target other versions of Windows.
#endif

////////////////////////////////////  For CRichToolTipCtrl
#include <afxdisp.h>
#include <afxole.h> // for COleClientItem
#define TTM_ADJUSTRECT          (WM_USER + 31)
#define TTP_STANDARD 1
#define TTSS_NORMAL 1
//////////////////////////////////////

#define WINVER 0x0600
#define _WIN32_WINNT 0x0600

#define VC_EXTRALEAN		
#include <afxwin.h>   
#include <afxext.h>   
#ifndef _AFX_NO_AFXCMN_SUPPORT
#include <afxcmn.h>			
#endif 
#include <winnls.h>	
#include <typeinfo.h> 
#include <math.h>
#include <string.h>
#include <omp.h>
#include <afxcolordialog.h>

#pragma warning( disable : 4244 4305 4239)  //conversion from 'double' to 'float', from 'const double' to 'const float'

// ========= STL
#include <map>
#include <list>
#include <vector> 
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <Windows.h>
#include <stdarg.h>
#include <algorithm>
#include <Shobjidl.h>

using namespace std;

//#include <afxdlgs.h>
//=== Подключение заголовков библиотек OpenGL
#include <gl/gl.h>
#include <gl/glu.h>
#include <GL/glaux.h>
//#include <GL/glut.h>

// podkluchenie vneshnix bibliotek
#pragma comment(lib,"OpenGL32.lib")
#pragma comment(lib,"Glu32.lib")
#pragma comment(lib,"GLAUX.LIB")
//#pragma comment(lib,"Glut32.lib")


typedef double DBL;    
//typedef float DBL; 
#include "atlimage.h" 
#include "resource.h"  

#include "C:\basic\_vt.h"
#include "C:\basic\_MyMFC.h"

//#include <cuda.h>
//#include <cuda_runtime_api.h>

//#include <mgl2/mgl.h>
#include <mgl/mgl_w.h>

//#include "myCuda.cuh"

#include "const.h"
#include "codes.h"
#include "matrix.h"
#include "helper.h"
#include "data.h"
#include "setup.h"
#include "filesystem.h"
#include "objects.h"
#include "problem.h"

#include "draw3.h"
#include "pde.h"
#include "pdesys.h"
#include "interactions.h"
#include "analysis3.h"
#include "draw.h"
#include "drawimg.h"
#include "move.h"

#include "save.h"

#include "dialog.h"
#include "Calculus.h"








