function [MC_wb, MC_db] = MC_conversion(m_water,m_solid)
n = length(m_water);
for i = 1:n
MC_wb(i) = m_water(i)/(m_solid(i)+m_water(i));
MC_db(i) = m_water(i)/m_solid(i);
end
end