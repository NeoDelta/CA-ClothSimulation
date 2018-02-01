#version 330 core
in vec3 Normal;
in vec3 Position;
out vec4 color;

uniform vec3 cameraPos;
uniform vec4 FragColor;

void main()
{             
    vec3 I = normalize(Position - cameraPos);
    color = FragColor;
}