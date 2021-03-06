function [ data, kr, kz, time ] = load_2D_data( filename, variablename )
%LOAD_2D_DATA load a 2D variable stored in a hdf5 result file from HeLaZ
    time     = h5read(filename,'/data/var2d/time');
    kr       = h5read(filename,['/data/var2d/',variablename,'/coordkr']);
    kz       = h5read(filename,['/data/var2d/',variablename,'/coordkz']);
    
    dt    = h5readatt(filename,'/data/input','dt');
    nsave = h5readatt(filename,'/data/input','nsave_5d'); 
    cstart = time(1)/dt/nsave;
    
    data     = zeros(numel(kr),numel(kz),numel(time));
    for it = 1:numel(time)
        tmp         = h5read(filename,['/data/var2d/',variablename,'/', num2str(it+cstart,'%06d')]);
        data(:,:,it) = tmp.real + 1i * tmp.imaginary;
    end

end

