function validate_and_save_response_objects(output_name, response_data)
%VALIDATE_AND_SAVE_RESPONSE_OBJECTS Save compact material response tensors.
%   The response-object convention is fixed to meV, nm, and nm^{-1}.

required_fields = { ...
    'material_name', ...
    'model_name', ...
    'version', ...
    'mu_grid', ...
    'S_photo', ...
    'Omega_z_photo', ...
    'D_photo', ...
    'J_photo', ...
    'units', ...
    'index_convention', ...
    'definitions', ...
    'parameters', ...
    'grid' ...
};

for i = 1:numel(required_fields)
    if ~isfield(response_data, required_fields{i})
        error('response_data is missing required field "%s".', required_fields{i});
    end
end

if ~strcmp(response_data.version, 'v1_response_objects')
    error('response_data.version must be "v1_response_objects".');
end

required_units.energy = 'meV';
required_units.length = 'nm';
required_units.momentum = 'nm^{-1}';
required_units.S = 'meV';
required_units.Omega_z = 'dimensionless';
required_units.D = 'nm';
required_units.J = 'meV nm';

unit_names = fieldnames(required_units);
for i = 1:numel(unit_names)
    name = unit_names{i};
    if ~isfield(response_data.units, name)
        error('response_data.units is missing required field "%s".', name);
    end
    if ~strcmp(response_data.units.(name), required_units.(name))
        error('response_data.units.%s must be "%s".', name, required_units.(name));
    end
end

Nmu = numel(response_data.mu_grid);
assert_size(response_data.S_photo, [Nmu 2 2], 'S_photo');
assert_size(response_data.Omega_z_photo, [Nmu 1], 'Omega_z_photo');
assert_size(response_data.D_photo, [Nmu 2], 'D_photo');
assert_size(response_data.J_photo, [Nmu 2 2 2], 'J_photo');

if ~isfield(response_data.definitions, 'D') || ...
        isempty(strfind(response_data.definitions.D, '- sum_n int_k Omega_n'))
    error('response_data.definitions.D must document D^i = - int Omega_xy partial_i f.');
end

save(output_name, '-struct', 'response_data');
end

function assert_size(value, expected, field_name)
actual = size(value);
if numel(actual) < numel(expected)
    actual(end+1:numel(expected)) = 1;
end
actual = actual(1:numel(expected));

if any(actual ~= expected)
    error('%s has size [%s], expected [%s].', ...
        field_name, num2str(actual), num2str(expected));
end
end
