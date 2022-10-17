using NCDatasets
using DataFrames

## call daily nutrient production data set, npp.
dp = NCDataset("/mnt/sda5-ext4/uci/oceanData/NPP_VAR_DATA/global-analysis-forecast-bio-001-028_lev1_8.nc");
dp2 = NCDataset("/mnt/sda5-ext4/uci/oceanData/NPP_VAR_DATA/global-analysis-forecast-bio-001-028_lev9_16.nc");
dp3 = NCDataset("/mnt/sda5-ext4/uci/oceanData/NPP_VAR_DATA/global-analysis-forecast-bio-001-028_lev17_24.nc");
dp4 = NCDataset("/mnt/sda5-ext4/uci/oceanData/NPP_VAR_DATA/global-analysis-forecast-bio-001-028_lev25_32.nc");
#[dp1["latitude"][:], dp2["latitude"][:], dp3["latitude"][:], dp4["latitude"][:]]
#[dp1["longitude"][:], dp2["longitude"][:], dp3["longitude"][:], dp4["longitude"][:]]
#[dp1["depth"][:], dp2["depth"][:], dp3["depth"][:], dp4["depth"][:]]
#[dp1["time"][:], dp2["time"][:], dp3["time"][:], dp4["time"][:]]
# latitude is from 35.0 -to- 45.0 with the stepsize of 0.25
lat_npp = dp["latitude"][:]
ny_npp  = length(lat_npp)
# longitude -80.0 -to- -5.0 with the stepsize of 0.25
lon_npp = dp["longitude"][:]
nx_npp  = length(lon_npp)
# depth from 0.494025 -to- 5792. with non-uniform vertical levels
depth_01_08=dp["depth"][:]
depth_09_16=dp2["depth"][:]
depth_17_24=dp3["depth"][:]
depth_25_32=dp4["depth"][:]
depth_npp = vcat(depth_01_08,depth_09_16,depth_17_24,depth_25_32)
nz_npp  = length(depth_npp)
depth_01_08=nothing;
depth_09_16=nothing;
depth_17_24=nothing;
depth_25_32=nothing;
# time is from 2019-01-01 -to- 2020-06-30 with the stepsize of 01 days
time_npp = dp["time"][:]
nt_npp  = length(time_npp)
# nppv = total Primary Production of Phyto (net primary production of biomass
# expressed as carbon per unit volume in sea water)
# data dimensions: longitude × latitude × depth × time = 301 × 41 × 11 × 547
# nppv = dp["nppv"][:,:,:,:]
nppv = Array{Union{Missing, Float32}}(missing,nx_npp,ny_npp,nz_npp,nt_npp)
nppv[:,:,1:8,:]   = dp["nppv"][:,:,:,:]
nppv[:,:,9:16,:]  = dp2["nppv"][:,:,:,:]
nppv[:,:,17:24,:] = dp3["nppv"][:,:,:,:]
nppv[:,:,25:32,:] = dp4["nppv"][:,:,:,:]

close(dp),close(dp2),close(dp3),close(dp4)
## call daily mean total surface current data set, tc. Note that the velocities
# are give at z=0 m and z=15 m
# absolute geostrophic velocity + Ekman velocity
ds = NCDataset("/mnt/sda5-ext4/uci/oceanData/NPP_VAR_DATA/dataset-uv-nrt-daily_1660575193826.nc");
#ds = NCDataset("/mnt/sda5-ext4/uci/oceanData/NPP_VAR_DATA/dataset-uv-nrt-daily_1652713133894.nc");

# latitude is from 35.125 -to- 45.125 with the stepsize of 0.25
lat_tc = ds["latitude"][:]
ny_tc = length(lat_tc) # no. of grid points in latitude direction
# longitude is form -79.875 -to- -5.125 with the stepsize of 0.25
lon_tc = ds["longitude"][:]
nx_tc = length(lon_tc) # no. of grid points in longitude direction
# depth at z = 0m (surface) and z = 15m
depth_tc = ds["depth"][:]
nz_tc = length(depth_tc) # no. of grid points in depth
# time is from 2019-01-01 -to- 2020-06-30 with the stepsize of 01 days
time_tc = ds["time"][:]
nt_tc = length(time_tc) # no. of time steps: per day
# zonal component: absolute geostrophic velocity + Ekman velocity
u = ds["uo"][:,:,:,:]
# meridian component: absolute geostrophic velocity + Ekman velocity
v = ds["vo"][:,:,:,:]

# close netcdf file
close(ds)
## get the common horizontal region and time in both the data sets of npp and tc
# commom longitude points
nx = min(nx_npp,nx_tc)
if nx_tc <= nx_npp
    lon=lon_npp[1:nx]
else
    lon=lon_tc[1:nx]
end
# commom latitude points
ny = min(ny_npp,ny_tc)
if ny_tc <= ny_npp
    lat=lat_npp[1:ny]
else
    lat=lat_tc[1:ny]
end
# common time points
nt = min(nt_npp,nt_tc)
if nt_tc <= nt_npp
    time=time_tc[1:nt]
else
    time=time_npp[1:nt]
end

## construct the grid (required for plotting)
# get the grid thickness in each direction
dx = (lon[2]-lon[1]) # grid thickness along longitude
dy = (lat[2]-lat[1]) # grid thickness along latitude
# grid points: generate the 2d matrix containing grid point
XG = repeat(lon',length(lat),1)
YG = repeat(lat,1, length(lon))
# compute grid thickness matrix
DX = dx*ones(ny,nx)
DY = dy*ones(ny,nx)
# grid centers: generate the 2d matrix containing grid centers
XC = XG .+ DX./2
YC = YG .+ DY./2

## reshape nppv and sum the npp variable over all vertical levels
# reshape the npp from shape (nx,ny,nz,nt) -to- (ny,nx,nz,nt)
nz = length(nppv[1,1,:,1]) # get the vertical length of nppv
npp = Array{Union{Missing, Float32}}(missing,ny,nx,nz,nt)
for j = 1:nz
    for i = 1:nt
        npp[:,:,j,i] = nppv[1:nx,1:ny,j,i]'
    end
end
# sum of nppv along depth (vertical levels)
npp_sum_over_depth=sum(npp, dims=3)
NPP = reshape(npp_sum_over_depth,(ny,nx,nt))

## reshape the velocities and compute the divergence
# reshape the velocity component from shape (nx,ny,nz,nt) to (ny,nx,nz,nt)
nz = length(u[1,1,:,1]) # get the vertical length of nppv
U = Array{Union{Missing, Float32}}(missing,ny,nx,nz,nt)
V = Array{Union{Missing, Float32}}(missing,ny,nx,nz,nt)
for j = 1:nz
    for i = 1:nt
        U[:,:,j,i] = u[:,:,j,i]'
        V[:,:,j,i] = v[:,:,j,i]'
    end
end
# compute the divervence of the velocity field
DUDX = Array{Union{Missing, Float32}}(missing,ny,nx-1,nz,nt)
DVDY = Array{Union{Missing, Float32}}(missing,ny-1,nx,nz,nt)
DIV = Array{Union{Missing, Float32}}(missing,ny-1,nx-1,nz,nt)
r = 6371.0;  # radius of the earth
for i in 1:nt
    for j in 1:nz
        DUDX[:,:,j,i] = diff(U[1:end,1:end,j,i],dims=2) ./ DX[1:end,1:end-1]
        DVDY[:,:,j,i] = diff(V[1:end,1:end,j,i].*(cos.(DY)),dims=1) ./ DY[1:end-1,1:end]
        DIV[:,:,j,i] = (1 ./ (r*cos.(DY[1:end-1,1:end-1]))).*(DUDX[1:end-1,1:end,j,i] + DVDY[1:end,1:end-1,j,i])
    end
end

## consider the subset of the input variables (this cell is for testing and deugging purpose)
mlat=1:ny-1;
mlon=1:nx-1;
# divergence at z = 0 meter
DIV_00 = DIV[mlat,mlon,1,:]
# divergence at z = 15 meter
DIV_15 = DIV[mlat,mlon,2,:]
# sum of net primary production (nppv) at z = 0.4 m -to- 11.81 m
N = NPP[mlat,mlon,:]

## start linear regression
# preallocate the matrix for R-squared values
R2 = zeros(length(mlat),length(mlon))
# adjust the time cycles for regression
nt = 2*365 # because data avalible day (replace this with nt to process at all avalible times in dataset)
T  = 365 # total number of days
t = 1:1:nt # time measure in days

for i in 1:length(mlat)
    for j in 1:length(mlon)
        # construct the predictive matrix A
        nh=3 # set no. of harmonics
        A = Array{Union{Missing, Float32}}(missing,nt,3+2*nh)
        A[:,1] = ones(nt,1) # column of ones
        A[:,2] = reshape(DIV_00[i,j,1:nt],nt,1) # divergence at z = 0 m
        A[:,3] = reshape(DIV_15[i,j,1:nt],nt,1) # divergence at z = 15 m
        l = 3
        for m in 1:nh
            A[:,l+1] = cos.(2*pi*m*t/T)
            A[:,l+2] = sin.(2*pi*m*t/T)
            l = l+2
        end
        # output data set of the observed npp
        Y  = reshape(N[i,j,1:nt],nt,1)
        # neglect the missing data values (i.e. land points will be ignored)
        idxY = any(ismissing, Y)
        idxA = any(ismissing, A)
        if  (idxA==true || idxY==true)
            # do nothing
        else
            b = A\Y # compute the regression co-efficients
            YModel = A*b#A[:,1:3]*b[1:3]
            e = sum((Y - YModel).*(Y - YModel))  # variance
            R2[i,j] = 1 - (e ./ sum(Y.*Y)) # compute r-squared
        end
    end
end

## plot R-squared
using Plots
gr()
xg = lon[mlon] # longitude grid centers
yg = lat[mlat] # latitude grid centers
minR2 = minimum(R2)
maxR2 = maximum(R2)
p1= contourf(xg,yg,R2,linewidth=0, levels=minR2:0.08:maxR2,showlabels = true,size=(1300,1000),xtickfont=18,ytickfont=18)
#p2=contourf(xg,yg,N[:,:,2],linewidth=0,levels=50,size=(1000,700),xtickfont=18,ytickfont=18)
#p1= contourf(xg,yg,R2)
plot(p1)

## animate npp
#xg = lon[mlon] # longitude grid centers
#yg = lat[mlat] # latitude grid centers
#begin
#    anim = @animate for i=1:100 # at time (In day)
#    #GR.setlinewidth(0)
#    contourf(xg,yg,N[:,:,i],linewidth=0,levels=50,size=(1000,700),xtickfont=18,
#    ytickfont=18)
#    #title!(string(Date(time[i])))
#    end
#    gif(anim, "/tmp/div.gif", fps = 1)
#end
