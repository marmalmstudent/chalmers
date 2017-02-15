clear
close all

%************** Parameters used in simulation (can be changed by the user)
N=256; % matrix size for the sampled fields in the optical system (Galilei�s telescope + eye)
lambda_simulation=500e-9; % simulation wavelength

fraction_ps_included=2/10000; % Maximally, each pixel in Saturn_HA3.jpg is one point source. But this would take far too long time
                              % for a normal computer. So this parameter tells how large fraction of the maximum number of point sources
                              % that is actually accounted for. (The pre-set value is too low as you will notice!)
D_sat=120e3*1e3; % diameter of Saturn (without rings)

f_obj=1.00; % focal length of objective lens
f_ocu=-5e-2; % focal length of ocular lens
f_eye=18e-3; % focal length of eye lens

L_sat=1400*1e6*1e3; % distance from Saturn to objective lens (Earth)
L_obj_ocu= 0.95; % distance from objective lens to ocular lens *CODE MISSING*
L_ocu_pupil=1e-2; % distance from ocular lens to pupil (and eye lens)
L_pupil_retina=f_eye; % distance from eye lens to retina (eye is assumed to be air-filled)

D_obj=2e-2; % diameter of objective lens
            % diameter of ocular lens is assumed large enough not to cut off any field
D_pupil_eye=3e-3; % diameter of eye pupil

size_of_numerical_window_retina=200e-6; % a "biologically" sensible default value if N=256 (the detector spacing in the retina is >= 1�m 
                                        % so there is not much point in having sampling distances much smaller)

k=2*pi/lambda_simulation;


%************** Reading in the intensity distribution in the "Saturn plane"
I_sat_raw=imread('Saturn_HA3.jpg');
I_sat=squeeze(I_sat_raw(:,:,2));
[number_rows,number_cols]=size(I_sat);
I_sat=double(I_sat).*double(I_sat>50); % A point source with too low intensity is set to zero (and will be omitted for speed)


%************** Initiate the "Saturn plane" (the sampling distance is set equal to the "true" size of one pixel in the image Saturn_HA3.jpg)
a_sat=D_sat/172; % Measuring in the image Saturn_HA3.jpg: the diameter of Saturn is approx 172 pixels
xvect_sat=-number_cols/2*a_sat:a_sat:(number_cols/2-1)*a_sat;
yvect_sat=-number_rows/2*a_sat:a_sat:(number_rows/2-1)*a_sat;
[xmat_sat,ymat_sat]=meshgrid(xvect_sat,yvect_sat);


%************** Initiate the plane of the objective lens
size_of_numerical_window_obj= 2e-2; % *CODE MISSING*
a_obj= size_of_numerical_window_obj/(N-1); % sampling distance *CODE MISSING*
xvect_obj=-N/2*a_obj:a_obj:(N/2-1)*a_obj;
yvect_obj=xvect_obj;
[xmat_obj,ymat_obj]=meshgrid(xvect_obj,yvect_obj);
rmat_obj=sqrt(xmat_obj.^2+ymat_obj.^2);
transmission_function_obj= exp(-1j*k*rmat_obj.^2/(2*f_obj)).*(rmat_obj < D_obj/2); % including cut-off of peripheral field *CODE MISSING*


%************** Initiate the plane of the ocular lens
size_of_numerical_window_ocu= 1.99e-3; % *CODE MISSING*
a_ocu= size_of_numerical_window_ocu/(N-1); % sampling distance *CODE MISSING*
xvect_ocu=-N/2*a_ocu:a_ocu:(N/2-1)*a_ocu;
yvect_ocu=xvect_ocu;
[xmat_ocu,ymat_ocu]=meshgrid(xvect_ocu,yvect_ocu);
rmat_ocu=sqrt(xmat_ocu.^2+ymat_ocu.^2);
transmission_function_ocu= exp(-1j*k*rmat_ocu.^2/(2*f_ocu)); % *CODE MISSING*


%************** Initiate the plane of the pupil/eye lens (this version is for the telescope simulation)
size_of_numerical_window_eyelens= 2e-3; % (this is for the telescope simulation, a 
                                    % different size of numerical window is suitable 
                                    % for the naked eye simulation) *CODE MISSING*
a_eyelens= size_of_numerical_window_eyelens/(N-1); % sampling distance *CODE MISSING*
xvect_eyelens=-N/2*a_eyelens:a_eyelens:(N/2-1)*a_eyelens;
yvect_eyelens=xvect_eyelens;
[xmat_eyelens,ymat_eyelens]=meshgrid(xvect_eyelens,yvect_eyelens);
rmat_eyelens=sqrt(xmat_eyelens.^2+ymat_eyelens.^2);
transmission_function_eyelens= exp(-1j*k*rmat_eyelens.^2/(2*f_eye)).*(rmat_obj < D_pupil_eye/2); % *CODE MISSING*


%************** Initiate the plane of the retina
a_retina= 2e-4/(N-1); % sampling distance *CODE MISSING*
xvect_retina=-N/2*a_retina:a_retina:(N/2-1)*a_retina;
yvect_retina=xvect_retina;
[xmat_retina,ymat_retina]=meshgrid(xvect_retina,yvect_retina);
rmat_retina=sqrt(xmat_retina.^2+ymat_retina.^2);


%************** Display the object, with your point source positions indicated
figure(1)
imagesc(xvect_sat/1e6/1e3,yvect_sat/1e6/1e3,I_sat)
colormap(gray)
axis('equal')
xlabel('milj km')
ylabel('milj km')
hold on
title(['Saturn with point sources indicated. Distance from Earth = ' num2str(round(L_sat/3e8/3600*10)/10) ' light-hours']) 

x_positions_sat=reshape(xmat_sat,1,number_rows*number_cols);
y_positions_sat=reshape(ymat_sat,1,number_rows*number_cols);
I_positions_sat=reshape(I_sat,1,number_rows*number_cols);

my_pixel_number=0;
for pixel_number=1:number_rows*number_cols
    if rem(pixel_number,round(1/fraction_ps_included))==0
       if I_positions_sat(pixel_number)>0
           my_pixel_number=my_pixel_number+1;
           my_x_pos_vect(my_pixel_number)=x_positions_sat(pixel_number);
           my_y_pos_vect(my_pixel_number)=y_positions_sat(pixel_number);
           my_I_sat_vect(my_pixel_number)=I_sat(pixel_number); % The intensity of the point source is assumed to be the brightness in the image pixel 
           plot(my_x_pos_vect(my_pixel_number)/1e6/1e3,my_y_pos_vect(my_pixel_number)/1e6/1e3,'.','MarkerEdgeColor',[0 0 0 ],'Markersize',6) % Indicate positions of the point sources on Saturn�s surface
       end
    end   
end
number_of_point_sources_in_simulation=my_pixel_number;



%************************************************************************************************************
%***************************************************************************** Start the iterative simulation
%************************************************************************************************************


I_retina_total=zeros(N,N); % accumulated intensity on retina

for point_source_number=1:number_of_point_sources_in_simulation % one point source at a time
    
    disp(['Point source number ' num2str(point_source_number) ' out of ' num2str(number_of_point_sources_in_simulation)])
    
    x_source=my_x_pos_vect(point_source_number); % x-position of point source on Saturn
    y_source=my_y_pos_vect(point_source_number); % y-position of point source on Saturn
    I_source=my_I_sat_vect(point_source_number); % intensity (brightness) of point source on Saturn
    
    figure(1)
    plot(x_source/1e6/1e3,y_source/1e6/1e3,'.','MarkerEdgeColor',[1 1 1 ]*0.95,'Markersize',7) % Indicate a source that turns on
    
        
    r_center= 1e-2; % distance from point source to the origin (0,0) in the plane of the objective lens *CODE MISSING*
    E_obj_in= ones(N,N); % incident field in (the sampling points in) the plane of the objective lens, from the point source *CODE MISSING*
    E_obj_out= E_obj_in.*transmission_function_obj; % field in the plane after the objective lens *CODE MISSING*
    I_obj_out=abs(E_obj_out).^2;
    figure(2)
    image(xvect_obj,yvect_obj,I_obj_out/max(max(I_obj_out))*64)
    colormap(gray)
    axis('equal')
    title('Intensity distribution after objective')
    
    
    [E_ocu_in,L1,L2,uvect,vvect,umat,vmat]= TSM(E_obj_out, a_obj, a_ocu, L_obj_ocu, lambda_simulation, 1); % field in the plane before the ocular lens *CODE MISSING*
    E_ocu_out= E_ocu_in.*transmission_function_ocu; % field in the plane after the ocular lens *CODE MISSING*
    I_ocu_out=abs(E_ocu_out).^2;
    figure(3)
    image(xvect_ocu,yvect_ocu,I_ocu_out/max(max(I_ocu_out))*64)
    colormap(gray)
    axis('equal')
    title('Intensity distribution after ocular')
    
    
    [E_eye_in,L1,L2,uvect,vvect,umat,vmat]= TSM(E_ocu_out, a_ocu, a_eyelens, L_ocu_pupil, lambda_simulation, 1); % field in the plane before the pupil/eye lens *CODE MISSING*
    E_eye_out= E_eye_in.*transmission_function_eyelens; % field in the plane after the pupil/eye lens *CODE MISSING*
    I_eye_out= abs(E_eye_out).^2;% Just to irritate you, you will need to fill it in yourself this time :) 
                % But seriously, this is the common way to display complex
                % fields, not the real part or anything like that! *CODE MISSING*
    figure(4)
    image(xvect_eyelens,yvect_eyelens,I_eye_out/max(max(I_eye_out))*64)
    colormap(gray)
    axis('equal')
    title('Intensity distribution after eye lens')
    
    
    [E_retina,L1,L2,uvect,vvect,umat,vmat]= TSM(E_eye_out, a_eyelens, a_retina, L_pupil_retina, lambda_simulation, 1.33); % field in the plane at the retina *CODE MISSING*
    I_retina=abs(E_retina).^2;
    figure(5)
    image(xvect_retina,yvect_retina,I_retina/max(max(I_retina))*64)
    colormap(gray)
    axis('equal')
    title('Intensity distribution at retina')
    
    
    I_retina_total=I_retina_total + I_retina; % the sum of intensity distributions from all point sources considered so far *CODE MISSING*
    figure(6)
    image(xvect_retina,yvect_retina,I_retina_total/max(max(I_retina_total))*64)
    colormap(gray)
    axis('equal')
    title('Intensity distribution at retina from all point sources considered so far')
    
    
    if point_source_number==1
        disp(['First point source completed. Place figures (1-6) so that they don�t overlap on screen. Then press any key to continue.'])
        pause
    end    
    pause(1)
end

