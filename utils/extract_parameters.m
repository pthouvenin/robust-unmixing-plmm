function [ Nsamples, Nlines, Nbands, interleave, offset, byte_order, data_type, wavelength, wavelength_unit ] = extract_parameters(PathName,FileName)
%Extract parameters from the header file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Error messages and warnings are generated when some of the keywors of interest are missing:

%Error : incomplete header file. No value found for the number of samples
%Error : incomplete header file. No value found for the number of lines.
%Error : incomplete header file. No value found for the number of bands.
%Error : incomplete header file. No interleave specification found.
%Error : incomplete header file. No byte order specified.
%Warning : no offset value specified. Default value assigned : offset = 0.
%Warning : no data type specified. Default value assigned : data type = double.
%Warning : no wavelength specified in the header file.

% Example of use :
% Filename = strcat(FileName, '.hdr'); %nom du fichier header dont les param�tres utiles sont extraits
% [Nsamples, Nlines, Nbands, interleave, offset, byte_order, data_type, wavelength, wavelength_unit] = extract_parameters(PathName,Filename);
%
% Irow = 1:Nlines;   % Indexes of rows
% Icol = 1:Nsamples; % Indexes of columns
% Iband = 1:Nbands;  % Indexes of bands
%
% FileName = strtok(FileName, '.'); %on enl�ve le .hdr de la cha�ne de caract�res initiale
% data_type = strcat(data_type,'=>','double'); %concat�nation des 2 strings
%
% data = multibandread(fullfile(PathName,FileName),[Nlines,Nsamples,Nbands],data_type,offset,interleave,byte_order,{'Column',Icol},{'Row',Irow},{'Band',Iband});

% Author: Pierre-Antoine Thouvenin, 2013.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Open / read .hdr file
fid = fopen(fullfile(PathName,FileName),'r');
X = fread(fid,[1,inf]);


% Look for appropriate keyword in the .hdr
id_samples = strfind(X,'samples');
id_lines = strfind(X,'lines');
id_bands = strfind(X,'bands');
id_wavelength = strfind(X,'wavelength');
id_interleave = strfind(X,'interleave');
id_offset = strfind(X,'offset');
id_byte_order = strfind(X,'byte order');
id_data_type = strfind(X,'data type');
id_wavelength_unit = strfind(X,'wavelength unit');


if isempty(id_samples)
    disp('Error : incomplete header file. No value found for the number of samples.');
    fclose(fid);
    return
else
    % first occurence of the word in the file
    pos_samples = id_samples(1,1);
    % extract related info
    status = fseek(fid,pos_samples,'bof'); % place cursor in the apporpriate position
    fseek(fid,9,'cof'); % skip word, spaces and = characters
    Nsamples = fscanf(fid,'%i'); % read data
end

if isempty(id_lines)
    disp('Error : incomplete header file. No value found for the number of lines.');
    fclose(fid);
    return
else
    pos_lines = id_lines(1,1);
    status = fseek(fid,pos_lines,'bof');
    fseek(fid,8,'cof');
    Nlines = fscanf(fid,'%i');
end

if isempty(id_bands)
    disp('Error : incomplete header file. No value found for the number of bands.');
    fclose(fid);
    return
else
    pos_bands = id_bands(1,1);
    status = fseek(fid,pos_bands,'bof');
    fseek(fid,8,'cof');
    Nbands = fscanf(fid,'%i');
end

if isempty(id_interleave)
    disp('Error : incomplete header file. No interleave specification found.');
    fclose(fid);
    return
else
    pos_interleave = id_interleave(1,1);
    status = fseek(fid,pos_interleave,'bof');
    fseek(fid,12,'cof');
    interleave = fscanf(fid,'%[bil bip bsq]'); % take only on of these 3 strings
end

if isempty(id_offset)
    disp('Warning : no offset value specified. Default value assigned : offset = 0.');
    offset = 0;
else
    pos_offset = id_offset(1,1);
    status = fseek(fid,pos_offset,'bof');
    fseek(fid,8,'cof');
    offset = fscanf(fid,'%i');
end

if isempty(id_data_type)
    disp('Warning : no data type specified. Default value assigned : data type = 5.');
    data_type = 5;
else
    pos_data_type = id_data_type(1,1);
    status = fseek(fid,pos_data_type,'bof');
    fseek(fid,11,'cof');
    data_type = fscanf(fid,'%i');
    switch data_type
        case 1
            data_type = 'int8';
        case 2
            data_type = 'int16';
        case 3
            data_type = 'int32';
        case 4
            data_type = 'float';
        case 5
            data_type = 'double';
        case 12
            data_type = 'uint16';
        otherwise
            data_type = 'float';
    end
end

if isempty(id_byte_order)
    disp('Error : incomplete header file. No byte order specified.');
    fclose(fid);
    return
else
    pos_byte_order = id_byte_order(1,1);
    status = fseek(fid,pos_byte_order,'bof');
    fseek(fid,12,'cof');
    byte_order = fscanf(fid,'%i');
    if byte_order == 0
        byte_order = 'ieee-le';
    else
        byte_order = 'ieee-be';
    end
end

if isempty(id_wavelength)
    disp('Warning : no wavelength specified in the header file.');
    wavelength=1:Nbands;
    wavelength=wavelength';
else
    [~,b]=size(id_wavelength); % if multiple occurences of he word "wavelength" in the file
    pos_wavelength = id_wavelength(1,b); % take last occurence of the word in the file
    status = fseek(fid,pos_wavelength,'bof');
    fseek(fid,13,'cof');
    
    % read wavelengths from the header file
    wavelength=zeros([Nbands 1]);
    for i = 1 : Nbands
        wavelength(i,1) = fscanf(fid,'%g');
        fseek(fid,2,'cof'); % skip separation comma
    end
end
%%
if isempty(id_wavelength_unit)
    disp('Warning : incomplete header file. No wavelength unit found.');
    wavelength_unit = 'NA';
    fclose(fid);
    return
else
    pos_wavelength_unit = id_wavelength_unit(1,1);
    status = fseek(fid,pos_wavelength_unit,'bof');
    fseek(fid,17,'cof');
    wavelength_unit = fscanf(fid,'%[Meters Millimeters Micrometers Nanometers]');
end

fclose(fid);

end
