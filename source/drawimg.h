#ifndef DRAWIMG_H
#define	DRAWIMG_H

using namespace Gdiplus;

class DrawImg
{
public:
	int width;
	VecL2 imgSize;
	VecF2 domSize;


public:
	DrawImg(Setup &s)
	{
		width = s.opt.imgW;
		domSize = s.spc.size;
		int sy = domSize.y / domSize.x * width;
		imgSize = VecL2(width, sy);
	}

public:
	void refresh(Setup &s)
	{
		width = s.opt.imgW;
		domSize = s.spc.size;
		int sy = domSize.y / domSize.x * width;
		imgSize = VecL2(width, sy);
	}

	void drawSimpleImage(string path, Setup &s, AtomList &al, PPL &plConns, PdeSys &ps)
	{
		Gdiplus::Bitmap *pDest = new Gdiplus::Bitmap(imgSize.x, imgSize.y, PixelFormat32bppRGB);
		Gdiplus::Graphics gr(pDest);

		Pen *pBlack10 = new Pen(Color::Black, 3);
		Pen *pRed1 = new Pen(Color::Red, 1);
		

		gr.Clear(Color::White);	

		if(s.opt.showFib && s.fib.enabled) drawFib(gr, s, ps);

		gr.DrawLine(pBlack10, 0,0,imgSize.x,0);
		gr.DrawLine(pBlack10, 0,imgSize.y,imgSize.x,imgSize.y);

		if(s.spc.prepArea) 
		{
			DBL x = getS(s.spc.sizePrepX);
			gr.DrawLine(pRed1, x,0,x,imgSize.y);
		}

		for(int i=0; i<al.obj.size(); i++)
			if(!al.obj[i].free) drawParticle(gr, s, al.obj[i]);

		drawConn(gr, s, plConns, al);

		CLSID pngClsid;
		GetEncoderClsid(L"image/png", &pngClsid);

		wstring ws(path.begin(), path.end());
		pDest->Save(ws.c_str(), &pngClsid, NULL);

		delete pDest;
		delete pBlack10, pRed1;
	}

private:
	int GetEncoderClsid(const WCHAR* format, CLSID* pClsid)
	{
		UINT  num = 0;		// number of image encoders
		UINT  size = 0;		// size of the image encoder array in bytes

		Gdiplus::ImageCodecInfo* pImageCodecInfo = NULL;

		Gdiplus::GetImageEncodersSize(&num, &size);
		if(size == 0)
			return -1;		// Failure

		pImageCodecInfo = (Gdiplus::ImageCodecInfo*)(malloc(size));
		if(pImageCodecInfo == NULL)
			return -1;		// Failure

		Gdiplus::GetImageEncoders(num, size, pImageCodecInfo);

		for(UINT j = 0; j < num; ++j)
		{
			if( wcscmp(pImageCodecInfo[j].MimeType, format) == 0 )
			{
				*pClsid = pImageCodecInfo[j].Clsid;
				free(pImageCodecInfo);
				return j;	// Success
			}    
		}

		free(pImageCodecInfo);
		return -1;  // Failure
	}

	inline VecL2 getP(VecF2 pos)
	{
		return VecL2(pos.x / domSize.x * imgSize.x, (1. - pos.y / domSize.y) * imgSize.y);
	}

	inline VecF2 getCoords(VecF2 i)
	{
		return VecF2(domSize.x * (i.x+0.5) / imgSize.x, (1. - (i.y+0.5) / imgSize.y) * domSize.y);
	}

	inline int getS(DBL h)
	{
		return h / domSize.x * imgSize.x;
	}

	void drawParticle(Graphics &gr, Setup &s, Atom &a)
	{
		Color col;

		if(a.isSPlt())
		{
			if(!s.opt.showPlt) return;

			if(a.oldPlt) col.SetFromCOLORREF(s.opt.cpla);
			else col.SetFromCOLORREF(s.opt.cplt);

			SolidBrush *b = new SolidBrush(col);

			VecL2 p = getP(a.pos - VecF2(s.atom.r, s.atom.r));
			DBL r = getS(s.atom.r * 2);
			RectF rect(p.x, p.y, r, -r); 
			gr.FillEllipse(b, rect);

			delete b;
		}
		else
		{
			if(!s.opt.showPla) return;

			col.SetFromCOLORREF(s.opt.cpla);
			SolidBrush *b = new SolidBrush(col);

			VecL2 p = getP(a.pos - VecF2(s.atom.rShow, s.atom.rShow));
			DBL r = getS(s.atom.rShow * 2);
			RectF rect(p.x, p.y, r, -r); 
			gr.FillEllipse(b, rect);

			delete b;
		}
	}

	void drawFib(Graphics &gr, Setup &s, PdeSys &ps)
	{
		Color col;
		col.SetFromCOLORREF(s.opt.cfib);
		SolidBrush *b = new SolidBrush(col);
		Pen *p = new Pen(col, 1);

		/*
		for(int i=0; i<imgSize.x-1; i++)
		{
			bool open = false;
			int beg = 0, end =0;
			int j=0;
			while(j<imgSize.y-1)
			{
				if(getU(i,j,s,ps)) 
				{
					if(!open)
					{
						beg = j;
						end = j;
						open = true;
					}
					else end = j;
				}
				else
				{
					if(open)
					{
						open = false;
						gr.DrawLine(p, i, beg, i, end);
					}
				}
				j++;
			}
			if(open)
			{
				end = imgSize.y-1;
				gr.DrawLine(p, i, beg, i, end);
			}
		}*/

		for(int j=0; j<imgSize.y-1; j++)
		{
			bool open = false;
			int beg = 0, end =0;
			int i=0;
			while(i<imgSize.x-1)
			{
				if(getU(i,j,s,ps)) 
				{
					if(!open)
					{
						beg = i;
						end = i;
						open = true;
					}
					else end = i;
				}
				else
				{
					if(open)
					{
						open = false;
						gr.DrawLine(p, beg, j, end, j);
					}
				}
				i++;
			}
			if(open)
			{
				end = imgSize.x-1;
				gr.DrawLine(p, beg, j, end, j);
			}
		}

		delete b,p;
	}

	bool getU(int i, int j, Setup &s, PdeSys &ps)
	{
		VecF2 c = getCoords(VecF2(i,j));
		return ps.pdeFp.getU(c) > s.plt.critFib;
	}

	void drawConn(Graphics &gr, Setup &s, PPL &plConns, AtomList &al)
	{
		Color cf1, cf2, cf3;
		cf1.SetFromCOLORREF(s.opt.cf1);
		cf2.SetFromCOLORREF(s.opt.cf2);
		cf3.SetFromCOLORREF(s.opt.cfib);
		Pen *pF1 = new Pen(cf1, 3);
		Pen *pF2 = new Pen(cf2, 3);
		Pen *pF3 = new Pen(Color::Black, 3);
		Pen *pen;

		for(vector<PP>::iterator iter = plConns.l.begin(); iter != plConns.l.end(); iter++)
		{
			VecL2 p1 = getP(al.obj[iter->a].pos);
			VecL2 p2 = getP(al.obj[iter->b].pos);

			
			if(al.obj[iter->a].oldPlt && al.obj[iter->b].oldPlt) pen = pF3;
			else if(iter->t > s.plt.F2t) pen = pF2;
			else pen = pF1;

			gr.DrawLine(pen, p1.x, p1.y, p2.x, p2.y);
		}

		delete pF1, pF2, pF3;
	}

};

#endif