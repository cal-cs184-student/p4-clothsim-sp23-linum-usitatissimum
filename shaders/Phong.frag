#version 330

uniform vec4 u_color;
uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

in vec4 v_position;
in vec4 v_normal;
in vec2 v_uv;

out vec4 out_color;

void main() {
  // YOUR CODE HERE
  float kd = 0.8;
  float ks = 0.5;
  float p = 128;
  float ka = 1.0;
  vec3 Ia = vec3(0.1, 0.1, 0.1);

  
  vec3 normal = normalize(v_normal).xyz;
  vec3 light_direction = normalize(u_light_pos - v_position.xyz);
  vec3 cam_direction = normalize(u_cam_pos - v_position.xyz);
  
  vec3 half_cam_direction = normalize(cam_direction + light_direction);

  float dist_light = pow(length(u_light_pos - v_position.xyz), 2.0);
  vec3 diffuse = kd * u_light_intensity / dist_light * max(0.0, dot(normal, light_direction));
  vec3 specular = ks * u_light_intensity / dist_light * pow(max(0.0, dot(normal, half_cam_direction)), p);
  out_color = vec4(ka * Ia + diffuse + specular, 1.0);
}

