 {~ Converts IES or EuLumDat photometry file to a table  ~}
 {~ or to the other format file as well ~}

program IES2tab;
   {2001 Jan Hollan, N.Copernicus Observatory and Planetarium }
(*  Copyright (C) 2001 Jan Hollan

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*)

uses dos, str_num, params, solar_ut, angles_o;
const
 InpFi:string[80]='';
 OutFIES:string[80]='';
 AngInRad:boolean=false;
 Heading:boolean=false;
 AsymY:boolean=false;
 AsymX:boolean=false;
 Flux  :real=0;
 FluxGlare:real=0;
 Flux80:real=0;
 Flux90:real=0;
 FluxU :real=0;
 FluxUq:array[1..4] of real = (0,0,0,0);
 FluxDown :real=0;
 FluxDownLow :real=0;
 IllumMax:real=-1;
 IllumMin:real=1E18;
 MaxLumInt80:real=-99;
 MaxLumInt90:real=-99;
 MaxLumIntM :real=-99;
 LineOnly:boolean=false;
 JustSegment:boolean=true;
 CalcUseful:boolean=false;
 SegmentSize:real=360;
 QMM:char=' ';
 QM8:char=' ';
 QM9:char=' ';
 xl:real=-0.5;
 xh:real= 1.5;
 yl:real=-3;
 yh:real= 3;
 CutoffType:array [0..7] of string[18] =
 ('Full_CutOff    ',
  'Fully_Shielded ',
  'CIE_CutOff     ',
  'IES_CutOff     ',
  'CIE_Semi-CutOff',
  'IES_Semi-CutOff',
  'Non-CutOff     ',
  'Unknown        ');
 CType:array [0..7] of string[3] =
 ('FCO','FS ','CO ','iCO',
  'SCO','iSC','NCO','unk');
 Select_Cutoff:boolean=false;
 SC:byte=1;
 OutFi:string[12]='sc.lis';
    VertAngMax=73; HoriAngMax=100;
    {VertAngMax=90; HoriAngMax=120;}
 cSky:array[4..6] of real=
 (10,16,24);
 dSky:array[4..6] of real=
 (-3,-3,-2.8);
 eSky:array[4..6] of real=
 (0.45,0.30,0.15);
 ZenExt:real=0.30;
  {default extinction in zenith is taken as 0.30 mag}
 Sky:byte=0;
 Albedo:real=0.10;
 LambDown:real=0.1108;
 LambDownLow:real=0.0528;
 LambDownSky:array[0..6] of real = (0.1108,
                                    0.1,0.1,0.1,
				    0.1475,
				    0.1303,
				    0.1093);
 LambDownSkyLow:array[0..6] of real = (0.0528,
                                    0.1,0.1,0.1,
				    0.0639,
				    0.0592,
				    0.0539);
 DefaultSky:boolean=true;
  {that with the default ZenExt=0.30 mag}
 InDegree=10;
 IncludeComments:boolean=false;
 EuLumDat: boolean=false;
 ForceCandMultiplier: boolean=false; 
  {true for wrong ies data giving cd/klm but not a true candela multiplier}
 MakeIES: boolean=false;
 LowAngle : real=30; {degrees above horizon considered for distant pollution}
 Distant : boolean=false;
 Tilt : boolean=false;
 Tilt1 : real=0; {around x-axis} 
 Tilt2 : real=0; {around y-axis}
 SegmentsTilt: byte=1;
 MaF : M3 = ((1,0,0), (0,1,0), (0,0,1));
 CT : byte = 0;
 MakeELD: boolean=false;
  MANUFAC:string[78]='';
  TEST   :string[78]='';
  LUMINAIRE:string[78]='';
  LUMCAT:string[78]='';
  DATE:string[78]='';
  LAMP:string[24]='';
  OTHERtilt:string[6]='0';

var     mess : string;
    ILine : string;
    Inpt,Outt,OutF: text;
    UpLimit: integer;
    i,j,k,HoriAngCount,VertAngCount,SegmentsI,segment: byte;
    LumFluxDklm, AngUnit, HoriAngInt, FluxElement, FluxElementGlare,
     Now, Last, CosNow, CosLast, CosBefLast, x,y, Illum, LnM, AerosFr,
     FluxXInt, sf, cf, c369, PolesAt80, VertAn, HoriAn, FluxElementHalf : real;
    PrFiel: array [1..13] of real;
    VertCosInt80,VertCosInt90,VertCosInt90d,VertCosInt75: real;
    VertAng, VertCosInt, VertTan, HoriCos, HoriSin, VertAngDownR:
     array [1..255] of real;
    HoriAng:
     array [0..255] of real;
    LumInt : array [1..VertAngMax,1..HoriAngMax] of real;
              {limitation is due to TP6 memory limitations,
               in newer Pascals the numbers can be higher.
               Alternatively, the numbers might be used immediately,
               just line by line, but the code would have to be changed.
               Of course, using pointers is also possible.}
    TooLong:boolean;
    D2r : real;
    EuSym : byte; {kind of EuLumDat symmetry}
    FDir: DirStr; FName: NameStr; FExt: ExtStr;
    laLen,laWid,laHei,laHei0,laHei90,laHei180,laHei270:real;
    Ma : M3; Ve: V3;


procedure Help;
begin
 writeln(mess);
 writeln(
'Program ies2tab converts IES photometric data to a table to standard output,',cl,
' or to a line with summary data (including the cutoff category).',cl,

' EULumDat format works as well, conversion to IES one is possible.',cl2,

' Its parameters are:',cl,
'-a#: Albedo of the ground (default 0.10)',cl,
'-c : write heading of Columns (just for an -l option)',cl,
{'-d : output in degrees (this is default)',cl,}
'-dp : compute sky luminance increment in distance instead of the overall one',cl,
'      (valid for the -l option); -di or -la30 are synonyma',cl,
'-e[{i|m}[<name]] : EuLumDat format',cl,
'   -e as an input;',cl, 
'   -ei            makes an almost-ies file too (default name *.ies)',cl, 
'               (almost: just due to letting >132 characters on a line)',cl,
'   -em as an output (default name *.ldt)',cl, 
'-f : assume the data giving the Full space angle of outcoming light',cl,
'-h : Help',cl,
'-i#: Indicatrix type (default 0, P.Cinzano (2000), 4..6 CIE types also possible)',cl,
'-ic: Include Comment lines in the output table',cl,
'-l : one Line output only',cl,
'-la# : limiting angular height for computing the increment of sky luminance',cl,
'       in distance, default 30 (degrees);',cl,
'       this increment is then given for the -l option',cl,
'-m : compute candela multiplier (if it is erroneously set to 1)',cl,
'-r : output angles in Radians (to ease polar plots by gnuplot)',cl,
'-s[#:[<name>]] : write filenames of <= cutoff categories to <name>',cl,
'    default limit is ',SC,' (0=FCO, 1=FS, 2=CIE-CO, 4=IES-CO, 5=CIE-SCO...),',cl,
'    default name is ',OutFi,cl,
'-t[y]# : tilt of the luminaire / 1 degree, around y-axis (perp. to the road)',cl,
' tx# : its tilt around x-axis (may be OK for long sloped roads)',cl, 
'-u[#:#:#:#] : compute illumination of a rectangle of',cl,
'    xmin:xmax:ymin:ymax luminaire heights, zero being below the luminaire',cl,
'    x axis going to 0 degrees and y to 90 degrees',cl,
'    (default ',xl:4:2,':',xh:4:2,':',yl:4:2,':',yh:4:2,' pole heights);',cl,
'-z#: extinction of light in Zenith / 1 mag (default 0.30)',cl2,
'( (C) Jan Hollan, N.Copernicus Observatory and Planetarium in Brno, 2001;'+cl+
' subject to the GNU General Public License, http://www.gnu.org/copyleft;'+cl+
' source code available at http://astro.sci.muni.cz/pub/hollan/programmes)');
 halt;
end;

procedure process_parameters;
begin
 for i:=1 to ParamCount do
  begin
   Par:=ParamStr(i);
   if Par[1] in [{'/',}'-'] then
    if length(Par)>1 then
     case Par[2] of
      'a','A': Albedo:=ss2r(Par,3);
      'c','C': Heading:=true;
      'd','D': if length(Par)>2 then
                if Par[3] in ['p','P','i','I'] then 
		 Distant:=true
      		else
	       else 	
    		AngInRad:=false;
      'e','E': begin
                 EuLumDat:=true;
		 if length(Par)>2 then
                 if Par[3] in ['i','I'] then 
                  begin
		   MakeIES:=true; 
		   if length(Par)>3 then
		    OutFIES:=copy(Par,4,length(Par)-3)
		   else
		    OutFIES:='';
		  end
		 else  
                 if Par[3] in ['m','M'] then 
                  begin
                   EuLumDat:=false;
		   MakeELD:=true; 
		   if length(Par)>3 then
		    OutFIES:=copy(Par,4,length(Par)-3)
		   else
		    OutFIES:='';
		  end
	       end;
      'f','F': JustSegment:=false;
      'h','?': Help;
      '-'    : if length(Par)>2 then
                if Par[3]='h' then Help;
      'i','I': if length(Par)>2 then
               if Par[3] in ['c','C'] then IncludeComments:=true
	       else if Par[3]='-' then IncludeComments:=false
	       else
                begin
                 Sky:=ss2i(Par,3);
                 if not (Sky in [0,4..6]) then
                  begin
                   Mess:='Indicatrix type is to be one of 0, 4, 5, 6';
                   help;
                  end
		 else begin LambDown:=LambDownSky[Sky];
		            LambDownLow:=LambDownSkyLow[Sky];
			end;		    
                end;
      'l','L': if length(Par)>2 then
                if Par[3] in ['a','A'] then 
		 begin
		  LowAngle:=ss2r(Par,4);
		  DefaultSky:=false;
		  Distant:=true;
		 end 
                else 
	       else LineOnly:=true;
      'm','M': ForceCandMultiplier:=true;
      'r','R': AngInRad:=true;
      's','S': begin
                Select_Cutoff:=true;
		if length(Par)>2 then
		SC:=I_S(ItemStrD(':',1,copy(Par,3,length(Par))));
		if ItemCountD(':',Par)=2 then
		 OutFi:=ItemStrD(':',2,Par);
               end;
      't','T': if length(Par)>2 then
                begin 
                 Tilt:=true; 
		 if not (Par[3] in ['x','X','y','Y']) then 
		  begin
		   Tilt2:=ss2r(Par,3);
		   mat(Ma,2,Tilt2*D2r);
		   matl_mat(Ma,MaF);
{		   writeln(MaF[1,1]MaF[2,1],MaF[3,1],cl,
		           MaF[1,2],MaF[2,2],MaF[3,2],cl, 
			   MaF[1,3],MaF[2,3],MaF[3,3]);
}			   
		  end
		 else
		  if length(Par)>3 then
                   if Par[3] in ['y','Y'] then
		    begin
                     Tilt2:=ss2r(Par,4);
                     mat(Ma,2,Tilt2*D2r);
		     matl_mat(Ma,MaF);
		    end 
		   else    
                    begin
                     Tilt1:=ss2r(Par,4);
                     mat(Ma,1,Tilt1*D2r);
		     matl_mat(Ma,MaF);
		    end; 
         	end;
      'u','U': begin
                CalcUseful:=true;
                if length(Par)>2 then
		 if length(Par)=3 then 
		  if Par[3]='-' then CalcUseful:=false else
		 else
		 if ItemCountD(':',Par)=4 then
		  begin
                   xl:=R_S(ItemStrD(':',1,copy(Par,3,length(Par))));
		   xh:=R_S(ItemStrD(':',2,Par));
		   yl:=R_S(ItemStrD(':',3,Par));
		   yh:=R_S(ItemStrD(':',4,Par));
		  end;
	       end;
      'z','Z': begin
                ZenExt:=ss2r(Par,3); DefaultSky:=false;
                if ZenExt<0.15 then
                 begin
                  Mess:='Zenith extinction is to be above 0.15 mag';
                  help;
                 end;
               end;
     end
    else
   else
    InpFi:=ParamStr(i);
  end;
end; {of process_parameters}

procedure AsResT(var f:text;s:string);
begin
 {$I-}
 assign(f,s);
 reset(f);
 if ioresult<>0 then
 begin
   mess:='There is no file '+s+' type the file name correctly!';
   Help;
 end;
 {$I+}
end;

procedure Quadrant(q:byte);
begin
 FluxUq[q]:=FluxUq[q] + FluxElement;
 if IllumMax<Illum then IllumMax:=Illum;
 if IllumMin>Illum then IllumMin:=Illum;
end;

procedure Useful(Hc,Hs:real);
      begin
       x:=VertTan[j]*Hc;
       y:=VertTan[j]*Hs;
       Illum:=sqrt(1+sqr(VertTan[j]));
       Illum:=LumInt[i,j]/(sqr(Illum)*Illum);
       if (xl<=x) and (x<=xh) and (yl<=y) and (y<=yh) then
        Quadrant(1);
       if AsymY and (SegmentsI>1) then
        begin
         if (xl<=x) and (x<=xh) and (-yh<=y) and (y<=-yl) then
          Quadrant(4);
	 if (SegmentsI=4) and AsymX then
	  if (-xh<=x) and (x<=-xl) and (-yh<=y) and (y<=-yl) then
	   Quadrant(3);
	end;
       if AsymX and (SegmentsI=4) then
        begin {to be added: reversed l and h}
         if (-xh<=x) and (x<=-xl) and (yl<=y) and (y<=yh) then
          Quadrant(2);
       	end;
      end;



procedure FluxDownVertAng(f:real; var FluxX:real);
 {this procedure integrates the downward-dispersed proportion 
  of the dispersed light (for Sky=0, just its aerosol fraction,
  the gaseous fraction is always one half)
  it is no longer used as the approximating function 
  given in the next procedure FluxDownForm gives this part instead,
  as obtained by fitting the results of this procedure previously}
var i,i1,i2:integer; xd,x,cx,sx,y,part:real;
begin
 FluxX:=0;
 sf:=sin(D2r*f);
 i1:=10*InDegree;
 i2:=124*InDegree;
 for i:=1 to UpLimit do
  begin
   xd:=(1.0/InDegree)*i;
   x:=D2r*xd;
   cx:=cos(x);
   sx:=sin(x);
   y:=sqr(xd);
   if xd<=f+1E-4 then
 
    {all the cone surface of dispersed light goes upwards} 
    part:=1
   else
    if xd>=180-f-1E-4 then
     {all the cone surface of dispersed light goes downwards}
     part:=0
    else
     begin
 
      part:=sf*cx/sx;
      part:=0.5+arctan(part/sqrt(1-sqr(part)))/pi;
     end;
   if part<1 then
    begin
     if Sky=0 then
      if i<=i1 then
       y:=7.5*exp(-0.1249*y/(1+0.04996*y))
      else
       if i<=i2 then
        y:=1.88*exp(-0.07226*xd+0.0002406*y)
       else y:=0.025 + 0.015*sin(2.25*x-c369)
     else
      y:=1+cSky[Sky]*(exp(dSky[Sky]*x)-exp(dSky[Sky]*pi/2))
       + eSky[Sky]*sqr(cx);
     FluxX:=FluxX+y*sx*(1-part);
    end;
  end;
 FluxX:=2*pi*FluxX*D2r/InDegree/FluxXInt;
end;

procedure FluxDownForm(a:real; var FluxX:real);
begin
 case Sky of
  0: FluxX:=(0.08+0.4*exp(-a/24))*AerosFr + 0.5*(1-AerosFr);
  4: FluxX:= 0.342+0.158*exp(-a/27);
  5: FluxX:= 0.275+0.225*exp(-a/27);
  6: FluxX:= 0.184+0.316*exp(-a/29);
 end;
end;

procedure LambDownNorm;
var {i:integer;} j:byte; {xd,x,cx,sx,}y:real;
begin

 {indicatrix integral:}
 {Not computed any more, as the values are needed
  just in the FluxDownVertAng procedure using non-normalized indicatrices.
  The functions in FluxDownForm are normalized already. }
{
if Sky>0 then
begin
 Flux:=0;
 for i:=1 to UpLimit do
  begin
   xd:=(1.0/InDegree)*i;
   x:=D2r*xd;
   cx:=cos(x);
   sx:=sin(x);
   y:=1+cSky[Sky]*(exp(dSky[Sky]*x)-exp(dSky[Sky]*pi/2)) + eSky[Sky]*sqr(cx);
   Flux:=Flux+y*sx;
  end;
 FluxXInt:=2*pi*Flux*D2r/InDegree;
end
else
 FluxXInt:=1.0027;
}
{4: 19.6061, 5: 22.0707, 6: 26.5463 would be the results for Sky 4, 5 or 6}

LambDown:=0;
LambDownLow:=0;

for j:=0 to 90 do
begin
 y:=j;
 cf:=cos(D2r*y);
 sf:=sin(D2r*y);
 {FluxDownVertAng(y,Flux);}
 FluxDownForm(y,Flux); {the approximated curves, instead of computing them}
 h2air_mass(y);
 y:=(1-1/exp(LnM*ZenExt*Air_Mass)) {dispersed fraction}
       *sf*cf;
  {sf: projection of the size of the source, 
                 cf: length of the horizontal projection of the unit vector} 
 {if Sky=0 then
   Flux:= Flux*AerosFr + 0.5*(1-AerosFr);} {adding the gaseous dispersion
                              when using FluxDownVertAng instead of 
			        	 FluxDownForm only}
 y:=y*Flux;	
 LambDown:=LambDown + y;
 if j<=LowAngle then 
  LambDownLow:=LambDownLow + y;
end;
 LambDown:=2*LambDown*D2r;
 LambDownLow:=2*LambDownLow*D2r;
 
{writeln(LambDown:8:4);
}  {the line is here for testing purposes only}
{writeln(LambDownLow:8:4);}
  {the line is here for testing purposes only}

end;


procedure VertAngDown(b:byte;a:real);
 var DispFr, FluxD:real;
{ function f2b(r:real):real;
  begin
   f2b:=exp(r*LnM)
  end;
 } 
 begin
  a:=a-90;
  h2air_mass(a);
  DispFr:=(1-1/exp(LnM*ZenExt*Air_Mass));
{  if DefaultSky then
   if Sky=0 then
     VertAngDownR[b]:=DispFr*((0.08+0.4*exp(-a/24))*AerosFr + 0.5*(1-AerosFr))
   else
     VertAngDownR[b]:=DispFr*(0.275+0.225*exp(-a/27))
  else
   begin
   FluxDownVertAng(a,Flux);}
  FluxDownForm(a,FluxD);
{    if Sky=0 then
     VertAngDownR[b]:=DispFr*(Flux*AerosFr + 0.5*(1-AerosFr))
   else  }
     VertAngDownR[b]:=DispFr*FluxD
{   end  }
 end;

procedure FillELD;
 procedure keywcoded(s:string; var so:string);
  var p,ls:byte;
  begin
   p:=pos('['+s+']',Iline);
   ls:=length(s);
   if p>0 then
    so:=LTrim(copy(ILine,p+ls+2,length(ILine)-(p+ls+1)))
   else so:='';
  end;
begin
 if MakeELD then 
      begin
        keywcoded('MANUFAC',mess); if mess<>'' then MANUFAC:=mess;
        keywcoded('TEST',mess); if mess<>'' then TEST:=mess;
        keywcoded('LUMINAIRE',mess); if mess<>'' then LUMINAIRE:=mess;
        keywcoded('LUMCAT',mess); if mess<>'' then LUMCAT:=mess;
        keywcoded('DATE',mess); if mess<>'' then DATE:=mess;
        keywcoded('LAMP',mess); if mess<>'' then LAMP:=mess;
      end;
end;


procedure WriteEulumdat;
 const
  EqdV:boolean=true;
  EqdH:boolean=true;
  IntPianiC:byte=0;
  IntGamma:byte=0;
 var tomm,Hint,Vint:real; Mc: byte; 
  laLen,laWid: word; {redefined here to be integers}
   begin
     Mc:=HoriAngCount;    
     case SegmentsI of
      1: if HoriAngCount=1 then EuSym:=1 else EuSym:=0;
      2: begin EuSym:=3; Mc:=2*(Mc-1); end;
      4: begin EuSym:=4; Mc:=4*(Mc-1); end;
      else EuSym:=5; {incompatible with both true IES and ELD format, 
                      all angles should be computed and EuSym:=0 set}
       {I'm not sure which IES case corresponds to EuSym=2}
     end; 

{testing if planes are equidistant:}
Hint:=360.0/Mc;
 if not (HoriAng[1]=0) then EqdH:=false;
 if (HoriAngCount>1) and EqdH then 
  begin
   for j:=3 to HoriAngCount do 
    if round(10*(HoriAng[j]-HoriAng[j-1]))<> round(10*Hint) then EqdH:=false;
  end;

 if not (VertAng[1]=0) then EqdV:=false;
 if (VertAngCount>2) and EqdV then 
  begin
   Vint:=VertAng[2];
   for j:=3 to VertAngCount do 
    if round(10*(VertAng[j]-VertAng[j-1])) <> round(10*Vint) then EqdV:=false;
  end;

if EqdH and (10*round(Hint)=round(10*Hint)) and (EuSym<>1) then 
 IntPianiC:= round(Hint);
if EqdV and (10*round(Vint)=round(10*Vint)) then 
 IntGamma:= round(Vint);
{... and setting ev. the increments if yes}

     if round(PrFiel[7])=2 then tomm:=1e3 
     else tomm:=304.8;
     laLen:=round(abs(PrFiel[9]*tomm));
     laWid:=round(abs(PrFiel[8]*tomm));
     if laLen=0 then begin laLen:=laWid; laWid:=0; end;  
     
     write
     (OutF,
      MANUFAC,cl, {1}
      '3',cl,
 {2 type 1 2 3}
      EuSym,cl, {3 symmetry type}
      {HoriAngCount}Mc,cl, {4 horiz. angles}
      IntPianiC,cl,    {5 interval}
      VertAngCount,cl, {6 vert. angles}
      IntGamma,cl,    {7 interval}
      TEST,cl,   {8 measurement report number 78}
      LUMINAIRE,cl, {9 luminaire name 78}
      LUMCAT,cl, {10 luminaire number 78}
      InpFi,cl, {11 file name 8}
      DATE,cl, {12 date/user .le. 78 }
      '0',cl, {13 length or diam./mm 4 }
      '0',cl, {14 width/mm 4, 0 for circ.}
      '0',cl, {15 height/mm 4 }
      laLen,cl, {16 length or diam. of luminous area /mm 4 }
      laWid,cl, {17 width/mm 4, 0 for circ.}
      {laHei0}round(abs(PrFiel[10]*tomm)),cl, {18 height 0 /mm 4 }
      {laHei90}'0',cl, {19 height 90 /mm 4 }
      {laHei180}'0',cl, {20 height 180 /mm 4 }
      {laHei270}'0',cl, {21 height 270 /mm 4 }
      round(100*(Flux-Flux90)/Flux),cl, {22 down flux ratio /% 4 }
      Flux/LumFluxDklm/10:4:1,cl, {23 light out /% 4 }
      '1',cl, {24 conv. factor  6 }
      OTHERtilt,cl, {25 tilt  6 }
      '1',cl, {26 bulb(lamp) sets  4 }
       {writeln(i); sets  4 }
      round(PrFiel[1]),cl, {a lamps}
      LAMP,cl,  {b lamp type}
      round(PrFiel[2]*PrFiel[1]),cl, {c flux} 
      '0',cl,  {d temperature}
      '0',cl, {e rendering}
      round(PrFiel[13]),cl {f watts} 
     );
 
      for j:=1 to 10 do write(OutF,'0',cl); {DR}
      for i:=1 to HoriAngCount do
 write(OutF,round(HoriAng[i]),cl);

      if EuSym>0 then 
       case EuSym of
{        1: for i:=1 to 35 do
 write(OutF,'0',cl);
}
        2,3: for i:={1 to 17} HoriAngCount+1 to Mc do
 write(OutF,{'0'} {round(HoriAng[1+i-HoriAngCount])+180}
                  180+(i-HoriAngCount)*IntPianiC ,cl);
		  {EulumDat assumes equidistant planes, so the two
		   expressions should give identic result}
        4: for i:={1 to 26}   HoriAngCount+1 to Mc do
 write(OutF,{'0'} 90 +(i-HoriAngCount)*IntPianiC,cl);
       end; 	 

      for j:=1 to VertAngCount do write(OutF,VertAng[j]:5:1,cl);
   end;
 {of WriteEuLumDat}

procedure ReadEulumdat;
   begin
      if MakeIES then
       write(OutF,'IESNA:LM-63-1995',cl,
	'[MANUFAC]   ',Iline,cl);
      if IncludeComments then 	
       writeln(
	'# [MANUFAC]   ',Iline);
      readln(Inpt);
 {2 type 1 2 3}
      readln(Inpt,EuSym); {3 symmetry type}
      readln(Inpt,HoriAngCount {PrFiel[5]}); {4 horiz. angles};
       {writeln(HoriAngCount); horiz. angles}
      readln(Inpt); {5 interval}
      readln(Inpt,VertAngCount {PrFiel[4]}); {6 vert. angles}
       {writeln(VertAngCount); vert. angles}
      readln(Inpt); {7 interval}
      readln(Inpt,Iline); {8 measurement report number}
      if MakeIES then write(OutF,
	'[TEST]      ',Iline,cl);
      if IncludeComments then 	
       writeln(
	'# [TEST]      ',Iline);
      readln(Inpt,Iline); {9 luminaire name 78}
      if MakeIES then write(OutF,
        '[LUMINAIRE] ',Iline,cl);
      if IncludeComments then writeln(	
        '# [LUMINAIRE] ',Iline);
      readln(Inpt,Iline); {10 luminaire number 78}
      if MakeIES then write(OutF,
        '[LUMCAT]    ',Iline,cl);
     
      if IncludeComments then writeln(
        '# [LUMCAT]    ',Iline);
     	
      readln(Inpt); {11 file name 8}
      readln(Inpt,Iline); {12 date/user <=78 }
      if MakeIES then write(OutF,
        '[DATE]      ',Iline,cl);
     
      if IncludeComments then writeln(
        '# [DATE]      ',Iline);
     
      readln(Inpt); {13 length or diam./mm 4 }
      readln(Inpt); {14 width/mm 4, 0 for circ.}
      readln(Inpt); {15 height/mm 4 }
      readln(Inpt,laLen); {16 length or diam. of luminous area /mm 4 }
      readln(Inpt,laWid); {17 width/mm 4, 0 for circ.}
      readln(Inpt,laHei0); {18 height 0 /mm 4 }
      readln(Inpt,laHei90); {19 height 90 /mm 4 }
      readln(Inpt,laHei180); {20 height 180 /mm 4 }
      readln(Inpt,laHei270); {21 height 270 /mm 4 }
      if MakeIES then
       begin
        laHei:=laHei0;
	if laHei90>laHei then laHei:=laHei90;
	if laHei180>laHei then laHei:=laHei180;
	if laHei270>laHei then laHei:=laHei270;
       end;      	
      readln(Inpt,Iline); {22 down flux ratio /% 4 }
      if IncludeComments then 	
       writeln('# down given: ',Iline, '%');
      readln(Inpt,Iline); {23 light out /% 4 }
      if IncludeComments then 	
       writeln('# out given: ',Iline, '%');
      readln(Inpt); {24 conv. factor  6 }
      readln(Inpt,Iline); {25 tilt  6 }
        if MakeIES then write(OutF,
        '[OTHER] tilt: ',Iline,cl);
     
      if IncludeComments then writeln(
        '# tilt: ',Iline); 	
      readln(Inpt,i); {26 bulb(lamp) sets  4 }
       {writeln(i); sets  4 }
      if i>0 then
       begin
        for j:=1 to i do readln(Inpt,PrFiel[1]); {lamps}
        for j:=1 to i do readln(Inpt,Iline); {lamp type}
        if MakeIES then write(OutF,
        '[LAMP]      ',Iline,cl);
     
        if IncludeComments then writeln(
        '# [LAMP]      ',Iline);
     
        for j:=1 to i do readln(Inpt,PrFiel[2]); {flux} 
        for j:=1 to i do readln(Inpt); {temperature}
        for j:=1 to i do readln(Inpt); {rendering}
        for j:=1 to i do readln(Inpt,PrFiel[13]); {watts}
	PrFiel[3]:=PrFiel[2]/1000; {bulb flux / klm}
       end;
      for j:=1 to 10 do readln(Inpt); {DR}
  
      read(Inpt,HoriAng[1]);
      if HoriAngCount>1 then
       for i:=2 to HoriAngCount do
        begin
         readln(Inpt,HoriAng[i]);
         if HoriAng[i]<=HoriAng[i-1] then
          begin
           writeln('Not rising horizontal angles.'); halt;
          end;
        end;
      for j:=1 to VertAngCount do readln(Inpt,VertAng[j]);
      if EuSym>0 then 
       begin 
        HoriAngCount:=HoriAngCount div 2;
        if EuSym>3 then HoriAngCount:=HoriAngCount div 2;
	HoriAngCount:=HoriAngCount + 1;
	if EuSym=1 then HoriAngCount:=1;
       end; 	 
      if MakeIES then 
       begin 
        write(OutF,
	 '[OTHER]   Eulumdat file:',InpFi,cl,
         'TILT=NONE',cl,
         PrFiel[1]:1:0, PrFiel[2]/PrFiel[1]:7:0, PrFiel[3]:5:1, 
	 VertAngCount:3, HoriAngCount:3, ' 1 2 ');
	if laWid<1 then {circular}
	 write(OutF, 
	  -laLen/1e3:5:2, ' 0', laHei/1e3:5:2,cl)
        else
	 write(OutF,
	  laWid/1e3:5:2, laLen/1e3:5:2, laHei/1e3:5:2,cl);
	write(OutF,
	 '1 1 ',PrFiel[13]:4:0,cl);
	for i:=1 to VertAngCount do write(OutF, VertAng[i]:6:1); 
	 write(OutF,cl);
	for i:=1 to HoriAngCount do write(OutF, HoriAng[i]:4:0); 
         write(OutF,cl);
       end;	
   end;
 {of EuLumDat}


begin
 D2r:=pi/180;
 c369:=369*D2r;
 LnM:=Ln(10)/2.5;
 Mess:='';
 if ParamCount=0 then Help;

 process_parameters;


 if (MakeIES or MakeELD) and (OutFIES='') then
  begin
   FSplit(InpFi,Fdir,Fname,Fext);
   if MakeELD then
    OutFIES:=Fdir+Fname+'.ldt'
   else    
    OutFIES:=Fdir+Fname+'.ies';
   if FSearch(OutFIES,'')<>'' then 
    begin
     writeln(OutFIES,' exists!');
     MakeIES:=false;
     MakeELD:=false;
    end 
  end;
 if AngInRad then AngUnit:=180/pi else AngUnit:=1;
 AerosFr:=(ZenExt-0.1)/ZenExt;
 UpLimit:=InDegree*180-1;

 if CalcUseful then 
  begin
   PolesAt80:=sin(80*D2r)/cos(80*D2r);
   TooLong:=
       (abs(xl)>PolesAt80)   
    or (abs(xh)>PolesAt80)   
    or (abs(yl)>PolesAt80)   
    or (abs(xh)>PolesAt80);

  end;
 if LineOnly then
  if Heading then
   begin
    if not Distant then
     writeln(
     '#  % of Increase of Sky Luminance')
    else
     writeln(
     '#  % of Increase of Sky Luminance in Distant Places by light below',
       LowAngle:5:1,' degrees',cl,
     '#  (for the zenith luminance such an angle suits places up to',
       5/sin(LowAngle*D2r):4:0,' km distance)');
    writeln(
     '#   due to light going from the luminaire above horizon, as compared',cl,
     '#   with the luminance produced by the light dispersed from the ground',cl,
     '#   concerns the following situation:',cl,
     '#  Albedo =', Albedo:5:2,cl,
     '#  Zenith Extinction =', ZenExt:5:2,' mag',
        ' (i.e., direct light remaining', 100/exp(LnM*ZenExt):4:0,' %)',cl,
     '#  Indicatrix type =', Sky:1,' (0: acc. to P.Cinzano, 4..6: CIE sky types)',cl);
    if CalcUseful then
     begin
      writeln(
      '#  (Calculation of the useful fraction of the outcoming luminous flux',cl,
      '#   concerns rectangle of ',xl:4:2,':',xh:4:2,':',yl:4:2,':',yh:4:2,' pole heights)',cl);
      if TooLong then writeln(
      '#   -- but just within',PolesAt80:5:2,' pole heights',cl2)
      else writeln;
     end; 
    write(
     '#FiOut',
     ' 75-90',
     ' >= 80degrees ',
     '   >= 90degrees ',
     '     CutOff?  ',
     '      filename');
    if CalcUseful then
     writeln(
     '  useful/Out ',
     'Illum/lx(1m,1klm)')
    else
     writeln;
    write(
     '#  %  ',
     '  %  ',
     '  cd/klm',
     ' % Out',
     '  cd/klm',
     ' % Out');
    if not Distant then
     write(
     ' % I.L.P. ')
    else
     write(
     ' % I.D.L.P.');
    if CalcUseful then
     writeln(
     '                            %',
     '    max   av   min  av/min')
    else
     writeln;
    if InpFi='' then halt;
   end;

 AsResT(Inpt,InpFi);
 if MakeIES or MakeELD then 
  begin
   assign(OutF,OutFIES);
   rewrite(OutF);
  end; 

 readln(Inpt,ILine);
 if eof(Inpt) then begin close(Inpt); halt(20) end
  {served for DOS version in case of a LF-only ended lines in the ies file}
 else
  if not EuLumDat then
  begin
   if IncludeComments then writeln('# ',Iline);
     FillELD;
   if not (copy(ILine,1,5)='TILT=') then
    repeat
     readln(Inpt,ILine); FillELD;
     if IncludeComments then writeln('# ',Iline);
    until (copy(ILine,1,5)='TILT=') or eof(Inpt);
  end
  else {if EuLumDat} ReadEulumdat;

 if eof(Inpt) then
  begin
   writeln(InpFi,': No  TILT=  line, probably no ies file.');
   halt(1);
  end;

 if copy(ILine,1,12)='TILT=INCLUDE' then 
  for i:=1 to 4 do readln(Inpt); {
skipping 16 numbers}

if not EuLumDat then 
begin

 for i:=1 to 13 do read(Inpt,PrFiel[i]);
 if ForceCandMultiplier then
  begin
   LumFluxDklm:=1;
   PrFiel[3]:=PrFiel[1]*PrFiel[2]/1000;
  end
 else
  LumFluxDklm:=PrFiel[1]*PrFiel[2]/1000/PrFiel[3];
 VertAngCount:=round(PrFiel[4]);
 if VertAngCount<>PrFiel[4] then
  begin writeln('Non-integer vert. angles count:',PrFiel[4]:6:1); halt end
 else if VertAngCount<1 then
  begin writeln('Non-positive vert. angles count:',PrFiel[4]:6:1); halt end;
 HoriAngCount:=round(PrFiel[5]);
 if HoriAngCount<>PrFiel[5] then
  begin writeln('Non-integer hori. angles count:',PrFiel[5]:6:1); halt end
 else if HoriAngCount<1 then
  begin writeln('Non-positive hori. angles count:',PrFiel[5]:6:1); halt end;

end
else {if EuLumDat}
  LumFluxDklm:=1; 

 if HoriAngCount>HoriAngMax then
  begin writeln('Too many horizontal angles,',
   HoriAngCount:4,'>',HoriAngMax:4); halt end;
 if VertAngCount>VertAngMax then
  begin writeln('Too many vertical angles,',
   VertAngCount:4,'>',VertAngMax); halt end;

if not EuLumDat then
begin
 
 if IncludeComments then
  begin
   writeln('# Bulbs:', PrFiel[1]:2:0, ',', PrFiel[2]:7:0, ' lm each');
   if PrFiel[3]<>1 then 
    writeln('#  (luminous intensites divided by', PrFiel[3]:6:2,
                ' in the original ies file)');
   writeln('# Number of measured angles:', PrFiel[4]:4:0,' vertical,', 
                                           PrFiel[5]:4:0,' horizontal ones');
   writeln('# Photometric type:', PrFiel[6]:2:0,
            ',  Units type:', PrFiel[7]:2:0);
   if PrFiel[8]<0 then	    
    writeln('# Luminous openig diameter (and height?):',
     -1*PrFiel[8]:8:3,'  (', PrFiel[10]:8:3,' )')
   else
   writeln('# Luminous openig width, length and height:',
    PrFiel[8]:8:3, PrFiel[9]:8:3, PrFiel[10]:8:3);
   writeln('# Ballast factor:', PrFiel[11]:2:0,
            ',  Input power:', PrFiel[13]:5:0,' W');
   writeln(cl);
  end;

 for j:=1 to VertAngCount do read(Inpt,VertAng[j]);
 if PrFiel[6]>=2 then {changing latitudes to polar distances}
  begin
   for j:=1 to VertAngCount do VertAng[j]:=VertAng[j]+90;
   if not Tilt then 
    begin
     CalcUseful:=false;
     CT:=7;
    end;
  end;

end;

 if not DefaultSky then  {downward fraction of the lambertian upward emission}
  
  LambDownNorm;
          {is to be computed}   

 Last:=VertAng[1];
 CosLast:=cos(Last*D2r);
 CosBefLast:=1;
 for j:=2 to VertAngCount-1 do
  begin
   if VertAng[j]<=VertAng[j-1] then
    begin
     writeln('Not rising vertical angles.'); halt;
    end;
   Now:=VertAng[j-1]+(VertAng[j]-VertAng[j-1])/2;
   CosNow:=cos(Now*D2r);
   VertCosInt[j-1]:=CosLast-CosNow;
   if VertAng[j-1]=90 then
    begin
     if CosBefLast=1 then CosBefLast:=-CosNow;
     VertCosInt90d:=CosBefLast;
 {down to the glare zone}
     VertCosInt90:=-CosNow;
     VertAngDown(j-1,(Now+90)/2)
    end
   else
    if VertAng[j-1]=80 then
     VertCosInt80:=cos(80*D2r)-CosNow
    else  
     if VertAng[j-1]=75 then
      VertCosInt75:=cos(75*D2r)-CosNow
     else if VertAng[j-1]>90 then
      VertAngDown(j-1,VertAng[j-1]);
   Last:=Now;
   CosBefLast:=CosLast;
   CosLast:=CosNow;
  end;
 j:=VertAngCount;
 Now:=VertAng[j];
 CosNow:=cos(Now*D2r);
 VertCosInt[j]:=CosLast-CosNow;
 if VertAng[j-1]>90 then
  VertAngDown(j-1,VertAng[j-1]);
 if VertAng[j]>90 then
  VertAngDown(j,VertAng[j]);

if not EuLumDat then
begin

 read(Inpt,HoriAng[1]);
 if HoriAngCount>1 then
  for i:=2 to HoriAngCount do
   begin
    read(Inpt,HoriAng[i]);
    if HoriAng[i]<=HoriAng[i-1] then
     begin
      writeln('Not rising horizontal angles.'); halt;
     end;
   end;

end;

 if CalcUseful then
  begin
   for j:=1 to VertAngCount do
    if VertAng[j]<89 then
     VertTan[j]:=sin(VertAng[j]*D2r)/cos(VertAng[j]*D2r);
   if HoriAngCount=1 then
    for j:=1 to 180 do
     begin
     HoriSin[j]:=sin((j-1)*2*D2r);
     HoriCos[j]:=cos((j-1)*2*D2r);
    end
   else
    for j:=1 to HoriAngCount do
     begin
      HoriSin[j]:=sin(HoriAng[j]*D2r);
      HoriCos[j]:=cos(HoriAng[j]*D2r);
     end;
   AsymY:=abs(yl+yh)>0.01;
   AsymX:=abs(xl+xh)>0.01;
  end;

 {$I-}
 for i:=1 to HoriAngCount do
 begin
  for j:=1 to VertAngCount do
   begin
    read(Inpt,LumInt[i,j]);
    if ioresult=106 then begin close(inpt); halt(106); end;
     {served for the case of commas-ended numbers in the ies file}
    if MakeIES then write(OutF,LumInt[i,j]:4:0); 
   end;
  if MakeIES then write(OutF,cl);
 end; 
 if MakeIES then close(OutF);
 close(inpt);
 {$I+}

{Flux integration and luminous intensity maxima:}
 if HoriAngCount>1 then
  begin
   HoriAng[HoriAngCount+1]:=HoriAng[HoriAngCount]+
    (HoriAng[HoriAngCount]-HoriAng[HoriAngCount-1]);
   HoriAng[0]:=HoriAng[1]-
    (HoriAng[2]-HoriAng[1])
  end
 else
  begin
   HoriAng[2]:=HoriAng[1]+180;
   HoriAng[0]:=HoriAng[1]-180;
  end;

 if HoriAng[HoriAngCount+1]-HoriAng[0]>=359.9 then
  begin
   SegmentSize:=360;
   JustSegment:=false;
  end;
 if JustSegment then
  begin
   SegmentSize:=HoriAng[HoriAngCount]-HoriAng[1];
   if PrFiel[6]=2 then {if (SegmentSize<90+1E-4) and (SegmentSize>90-1E-4)}
    if HoriAng[1]=0 then SegmentSize:=180
    else SegmentSize:=360;
  end;
 if SegmentSize<89.9 then CalcUseful:=false;
 SegmentsI:=round(360/SegmentSize);
 if abs(360/SegmentSize - SegmentsI) >0.1 then
  begin
   JustSegment:=false;
   SegmentSize:=360;
   SegmentsI:=1;
  end;
 if (SegmentsI>1) and (not(SegmentsI in [2,4])) then
  CalcUseful:=false;
 if VertAng[1]>80 then CalcUseful:=false; {for uplights}

 if Tilt then
  begin
{   mat(Ma,2,Tilt2*D2r);}
   SegmentsTilt:=SegmentsI;
   if HoriAngCount=1 then
    begin
     SegmentsTilt:=36;
     SegmentSize:=10;
    end; 
  end;

 Last:=HoriAng[0]+(HoriAng[1]-HoriAng[0])/2;
 for i:=1 to HoriAngCount do
  begin
   Now:=HoriAng[i]+(HoriAng[i+1]-HoriAng[i])/2;
   HoriAngInt:=Now-Last;
   if HoriAngCount=1 then 
    if Tilt then  
     HoriAngInt:=10
    else 
     HoriAngInt:=360;
   Last:=Now;
   if JustSegment and ((i=1) or (i=HoriAngCount)) then
    HoriAngInt:=HoriAngInt / 2;
   for segment:=1 to SegmentsTilt do
   for j:=1 to VertAngCount do 
    begin
     if not Tilt then
      begin
       VertAn:=VertAng[j];
      end
     else
      begin
       {rotation, bringing the VertAn to latitude temporarily:}
        VertAn:=(VertAng[j]-90)*D2r;
	if Tilt and (SegmentsI<6) then  
	 HoriAn:=(HoriAng[i]*(2*(segment mod 2)-1)
	          + 180 * (segment div 3)
		  )*D2r
	else {no plane symmetry, just repeating segments}
	 HoriAn:=(HoriAng[i]+SegmentSize*(segment-1))*D2r;
{writeln(HoriAn/D2r:5:1, VertAn/D2r:5:1);    }
        P3s_Vector(HoriAn,VertAn,Ve); {V[3] to max F, V[1] to 0 L}
        mat_vect(MaF,Ve);
          {vector:= matrix * vector}
        Vector_P3s(Ve,HoriAn,VertAn);
        if CalcUseful then
	 begin 
	  HoriSin[j]:=sin(HoriAn);
          HoriCos[j]:=cos(HoriAn);
	 end; 
	HoriAn:=HoriAn/D2r;
	VertAn:=90+VertAn/D2r;
{writeln(HoriAn:5:1, VertAn:5:1);   }
      end;
     FluxElement:=LumInt[i,j] * HoriAngInt * VertCosInt[j];
     FluxElementHalf:=FluxElement/2;
     Flux:=Flux + FluxElement;
     if VertAn>=75 then
      if VertAn<=90 then
        begin {glare part}
         FluxElementGlare:=FluxElement;
         if VertAn=75 then {taking just the part above 75 degrees}
          if Tilt then
	   FluxElementGlare:=FluxElementHalf
	  else
	   FluxElementGlare:=LumInt[i,j] * HoriAngInt * VertCosInt75
	 else 
          if VertAn=90 then
           if Tilt then
            FluxElementGlare:=FluxElementHalf
	   else
            FluxElementGlare:=LumInt[i,j] * HoriAngInt * VertCosInt90d
	  else FluxElementGlare:=FluxElement; 
         if VertAn<=90 then
          FluxGlare:=FluxGlare + FluxElementGlare;
	end;  
     if VertAn>=80 then
      begin
       if VertAn=80 then {taking just the part above 80 degrees}
        if Tilt then
         FluxElement:=FluxElementHalf
        else
         FluxElement:=LumInt[i,j] * HoriAngInt * VertCosInt80;
       Flux80:=Flux80 + FluxElement;
       if LumInt[i,j]>MaxLumInt80 then MaxLumInt80:=LumInt[i,j];
       if VertAn>=90 then
        begin
	 if VertAn=90 then
          if Tilt then
           FluxElement:=FluxElementHalf
          else
           FluxElement:=LumInt[i,j] * HoriAngInt * VertCosInt90;
         Flux90:=Flux90 + FluxElement;
	 if Tilt then 
	  if VertAn=90 then
           begin
     	    if j<VertAngCount then
             Now:=(VertAng[j+1]-VertAng[j])/2
	    else    
             Now:=(VertAng[j]-VertAng[j-1])/2;
            VertAngDown(j,90+Now/2);
	   end
	  else
	   VertAngDown(j,VertAn);
	 FluxElement:=FluxElement*VertAngDownR[j]; {becomes downflux instead}
	 FluxDown:=FluxDown + {FluxElement*VertAngDownR[j]}
	                       FluxElement;
	 if VertAn<=LowAngle+90 then		       
	  FluxDownLow:=FluxDownLow + FluxElement;
         if LumInt[i,j]>MaxLumInt90 then MaxLumInt90:=LumInt[i,j];
	 if j=VertAngCount then
          if LumInt[i,j]>MaxLumIntM then MaxLumIntM:=LumInt[i,j];
        end;
      end
     else if CalcUseful then {assuming no useful illum above 80 degrees}
      if HoriAngCount=1 then
       begin
        FluxElement:=FluxElement/180;
        for 
	k:=1 to 180 do Useful(HoriCos[k],HoriSin[k])
       end	
      else
       Useful(HoriCos[i],HoriSin[i]);
    end;  
  end;

 if SegmentsTilt>1 then SegmentSize:=360;
 Flux  :=Flux  *360*D2r/SegmentSize;
 Flux80:=Flux80*360*D2r/SegmentSize;
 Flux90:=Flux90*360*D2r/SegmentSize;
 FluxGlare:=FluxGlare*360*D2r/SegmentSize;
 FluxDown:=FluxDown*360*D2r/SegmentSize;
 FluxDownLow:=FluxDownLow*360*D2r/SegmentSize;
 MaxLumInt80:=MaxLumInt80/LumFluxDklm;
 MaxLumInt90:=MaxLumInt90/LumFluxDklm;
 if CalcUseful then
  begin
   FluxU :=FluxUq[1];
   case SegmentsI of
    2: if AsymY then
        FluxU :=FluxU+FluxUq[4]
       else FluxU:=2*FluxU;
    4: begin
        if AsymY then
         FluxU :=FluxU+FluxUq[4]
        else FluxU:=2*FluxU;
        if AsymX then
         FluxU :=FluxU+FluxUq[2]
        else FluxU:=2*FluxU;
        if AsymX and AsymY then
         FluxU :=FluxU+FluxUq[3]
       end;	
   end;
   FluxU :=FluxU*D2r;
  end;

 if CT<>7 then      	
 if (MaxLumInt90<0.5) and (MaxLumInt80<100) then CT:=0
 else
  if MaxLumInt90<0.5 then CT:=1
  else
   if (MaxLumInt90<=10) and (MaxLumInt80<=30) then CT:=2
   else
    if (MaxLumInt90<=25) and (MaxLumInt80<=100) then CT:=3
    else
     if (MaxLumInt90<=50) and (MaxLumInt80<=100) then CT:=4
     else
      if (MaxLumInt90<=50) and (MaxLumInt80<=200) then CT:=5
      else CT:=6;

 if (VertAng[VertAngCount]<90) and (MaxLumInt80>0) then QM8:='?';
 if (VertAng[VertAngCount]=90) and (MaxLumInt90>0) then QM9:='?';
 if (VertAng[VertAngCount]<179.9) and (MaxLumIntM>0) then QMM:='?';

 if select_cutoff then
  if CT<=SC then
   begin
    {$I-}
     assign(outt,OutFi);
     append(outt);
     if ioresult<>0 then
     rewrite(outt);
    {$I+}
    write(outt,InpFi,' ');
    close(outt);
  end;

 if Flux/LumFluxDklm/10>100 then
  writeln(
 '# NONSENSE SUMMARY RESULTS of ies2tab: output flux is to be under 100 %!');

if MakeELD then
begin
 WriteEulumdat;
 for i:=1 to HoriAngCount do
  for j:=1 to VertAngCount do
    write(OutF,round(LumInt[i,j]),cl); 
 if MakeELD then close(OutF);
end;

 if LineOnly then
  begin
   write(
    Flux/LumFluxDklm/10:6:1,
    100*FluxGlare/Flux:5:1,
    MaxLumInt80:7:1,QM8,
    100*Flux80/Flux:5:1,
    MaxLumInt90:7:1,QM9,
    100*Flux90/Flux:5:1,QMM);
   if not Distant then 
    write(   
    100*FluxDown/((Flux-Flux90)*Albedo*LambDown):4:0,QMM)
   else
    write(
    100*FluxDownLow/((Flux-Flux90)*Albedo*LambDownLow):4:0,QMM);
   write(
    '  ',CutoffType[CT],' ',InpFi);
   if length(InpFi)<12 then
    for j:=1 to 12-length(InpFi) do write(' ');
   if CalcUseful then
    begin
     write(
      100*FluxU/Flux:7:1,
      IllumMax/LumFluxDklm:7:1,
      FluxU/LumFluxDklm/((xh-xl)*(yh-yl)):6:1,
      IllumMin/LumFluxDklm:5:1);
     if IllumMin/LumFluxDklm>0.1 then writeln(
       FluxU/((xh-xl)*(yh-yl))/IllumMin:6:1)
     else writeln
    end
   else
    writeln;
   halt
  end
 else
  write(
   '# Source file: ',InpFi,cl,
   '# Luminaire flux =',Flux:6:0);
  if PrFiel[3]<>1 then 
    writeln(' raw, for the given bulb(s) it would be', Flux*PrFiel[3]:6:0,' lx,')
  else writeln(' lx,');    
  writeln(    
   '#                 ',Flux/LumFluxDklm/10:6:1,' % of the bulb flux',cl,
   '# between 75 and 90: ',100*FluxGlare/Flux:4:1,' % of the luminaire flux',cl,
   '#  - this part causes just GLARE in case of road lighting and similar purposes',cl,                
   '# 80deg and above:  max',MaxLumInt80:6:1,QM8,'cd/1000lm ,',
      100*Flux80/Flux:5:1,' % of the luminaire flux', cl,
   '# 90deg and above:  max',MaxLumInt90:6:1,QM9,'cd/1000lm ,',
      100*Flux90/Flux:5:1,QMM,'% of the luminaire flux',cl,
   '# CutOff Type: ',CutoffType[CT],cl);
    writeln(
     '#  Increase of Sky Luminance due to light going',cl,
     '#    from the luminaire difectly above horizon, as compared with the',cl,
     '#    luminance produced by the light dispersed from the ground:',
       {100*FluxDown/((Flux-FluxDown)*Albedo*LambDown):4:0,QMM,'%',cl,}
                   {queer, Flux90 should be here and elsewhere 
		instead of FluxDown? For globes, there is a difference 
		-- this error was present before Oct 18, 2002}
     {'#    no, this is more accurate:',}
       100*FluxDown/((Flux-Flux90)*Albedo*LambDown):4:0,QMM,'%',cl,
     '#  Increase of Sky Luminance in Distant Places by light below',
     LowAngle:5:1,' degrees',cl,
     '#    due to light going from the luminaire directly above horizon:',
       100*FluxDownLow/((Flux-Flux90)*Albedo*LambDownLow):4:0,QMM,'%',cl,
     '#  (for the zenith luminance such an angle suits places up to',
       5/sin(LowAngle*D2r):4:0,' km distance)',cl,
     '#   The increases concern the following situation:',cl,
     '#    Albedo =', Albedo:5:2,cl,
     '#    Zenith Extinction =', ZenExt:5:2,' mag',
        ' (i.e., direct light remaining', 100/exp(LnM*ZenExt):4:0,' %)',cl,
     '#    Indicatrix type =', Sky:1,' (0: acc. to P.Cinzano, 4..6: CIE sky types)',cl,
     '#   (the downward-scattered part of lambertian uplight is',
           LambDown:7:4, ' then)',cl);
   if CalcUseful then
    begin
     write(
   '# Illumination of a rectangle of ',
    xl:3:1,':',xh:3:1,':',yl:3:1,':',yh:3:1,' pole heights');
    
 if TooLong then  writeln(
   '#   -- but just within',PolesAt80:5:2,' pole heights:')
     else writeln(':');
     write(
   '#  Per cent of luminaire output falling there:', 100*FluxU/Flux:7:1,cl,
   '#   for a unit case (glass at 1 m height, 1000 lm bulb) its illumination:',cl,
   '#  Maximum: ', IllumMax/LumFluxDklm:7:1,' lx',cl,
   '#  Average: ', FluxU/LumFluxDklm/((xh-xl)*(yh-yl)):6:1,' lx',cl,
   '#  Minimum: ', IllumMin/LumFluxDklm:6:1,' lx',cl);
     if IllumMin/LumFluxDklm>0.1 then writeln(
   '#  Average/Minimum: ',FluxU/((xh-xl)*(yh-yl))/IllumMin:6:1,cl)
     else writeln(cl);
    end;

 writeln(
'# The following table gives luminous intensities which would be produced',cl,
'# using a hypotetic bulb giving a luminous flux of 1000 lm (i.e., cd/klm):',cl);
  
 if HoriAngCount>1 then
  begin
   write('#  H:'); for i:=1 to HoriAngCount do write(HoriAng[i]:6:1); writeln;
   write('#V:  '); writeln;
  end;
 for j:=1 to VertAngCount do
  begin
   write(VertAng[j]/AngUnit:5:1);
   for i:=1 to HoriAngCount do write(LumInt[i,j]/LumFluxDklm:6:1);
   writeln;
  end;

end.
