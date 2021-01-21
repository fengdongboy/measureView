#ifdef GL_ES
// Set default precision to medium
precision highp int;
precision highp float;
#endif

     uniform  vec3 uLightPosition;
     uniform  vec3 uLightAmbient;
     uniform  vec3 uLightDiffuse;
     uniform  vec3 uLightSpecular;

     uniform  vec3 uMaterialAmbient;
     uniform  vec3 uMaterialDiffuse;
     uniform  vec3 uMaterialSpecular;
     uniform  float uShininess;
	 uniform  int renderType;

     varying  vec3 fNormal;
     varying  vec3 fPosition;

     void main()
     {
        vec3 light_direction = normalize(uLightPosition - fPosition);
        float diffuse_value = max(0.0,dot(light_direction,fNormal));
        vec3 reflect_direction = reflect(-light_direction,normalize(fNormal));
        vec3 eye_direction = normalize(-fPosition);
        float cosvalue_lighteye = max(0.0,dot(reflect_direction,eye_direction));
        float specular_value = 0.0;
        if(diffuse_value != 0.0)
                specular_value = pow(cosvalue_lighteye,uShininess);
        vec3 ambient_light = uLightAmbient * uMaterialAmbient;
        vec3 diffuse_light = diffuse_value * uLightDiffuse * uMaterialDiffuse;
        vec3 specular_light = specular_value * uLightSpecular * uMaterialSpecular;
        vec3 ads_light = ambient_light+diffuse_light+specular_light;
        vec3 final_color = min(ads_light,vec3(1.0));
        //gl_PointSize = 4.0;
		if (renderType == 1)
			 gl_FragColor = vec4(0.95,0.95,0,1.0);
		else if(renderType == 2)
			gl_FragColor = vec4(0.12,0.411,0.82,1.0);
		else if (renderType == 3)
			gl_FragColor = vec4(final_color,0.6);
		else if (renderType == 4)
		{
		    gl_FragColor = vec4(1.0,0,0,1.0);
		}
		else if (renderType == 5)
			gl_FragColor = vec4(0.0,1.0,0.0,1.0);
		else 
			gl_FragColor = vec4(final_color,1.0);
        //gl_FragColor=vec4(1.0,0.0,0.0,1.0);
     }
